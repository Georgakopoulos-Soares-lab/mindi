# Keep your GFF Clean & Tidy

__author__ = "Nikol Chantzi"
__email__ = "nmc6088@psu.edu"
__version__ = "1.0.1"

### > IMPORTS BEGIN 

import os
import tempfile
import pybedtools
from pybedtools import BedTool
from pathlib import Path
import subprocess
from typing import Optional, ClassVar
import pandas as pd

### < IMPORTS END

class GFFCleaner:
    
    GFF_FIELDS: ClassVar[list[str]] = [
                                       "seqID", 
                                       "source",
                                       "compartment",
                                       "start", 
                                       "end",
                                       "score",
                                       "strand",
                                       "phase",
                                       "attributes"
                                ]

    def __init__(self, 
                 tempdir: Optional[os.PathLike[str]] = None, 
                 all_compartments: Optional[list[tuple[str]]] = None,
                 bedtools_path: Optional[os.PathLike[str]] = None,
                 ) -> None:
        
        # CONFIGURATION

        if tempdir is None:
            tempdir = Path().cwd()
        else:
            tempdir = Path(tempdir).resolve()
            tempdir.mkdir(exist_ok=True)

        self.tempdir = tempdir
        pybedtools.set_tempdir(self.tempdir)
        pybedtools.set_bedtools_path(bedtools_path)
        
        if all_compartments is None:
            all_compartments = [
                                ("region", None),
                                 ("CDS", None),
                                 ("exon", None),
                                 ("gene", None),
                                 ("gene", "protein_coding"),
                                 ("gene", "non_coding"),
                            ]
        
        self.valid_compartments = {"region", "CDS", "exon", "gene"}
        self.all_compartments = all_compartments

    
    @staticmethod
    def parse_attributes(attributes: str) -> dict[str, str]:
        attributes = attributes.strip()
        if "=" not in attributes:
            return attributes

        splitted_attributes = attributes.split(";")
        parsed_attributes = dict()
        for attr in splitted_attributes:
            name, value = attr.split("=")
            parsed_attributes.update({name: value})

        return parsed_attributes

    @staticmethod
    def parse_biotype(attributes: str) -> str:
        if "gene_biotype=" not in attributes:
            return "."

        biotype = attributes.split("gene_biotype=")[1]

        if ";" in biotype:
            return biotype.split(";")[0]
        return biotype


    @staticmethod
    def reduce_biotype(biotype: str) -> str:
        if biotype == "protein_coding":
            return biotype
        elif biotype != "." and biotype is not None:
            return "non_coding"
        return "."

    
    def read_gff(self, gff: os.PathLike[str]) -> pd.DataFrame:
        gff_df = pd.read_table(
                            gff,
                            comment="#",
                            header=None, 
                            names=GFFCleaner.GFF_FIELDS,
                            dtype={
                                "start": int,
                                "end": int
                            }
                    )
        gff_df = gff_df[gff_df['compartment'].isin(self.valid_compartments)].reset_index(drop=True)
        return gff_df


    def read(self, gff: os.PathLike[str], add_exons: bool = True, reduce_biotype: bool = True) -> pd.DataFrame:
        gff = Path(gff).resolve()
        gff_name = gff.name.split('.gff')[0]
        
        unzipped_gff = None
        
        if add_exons:

            reading_temp_dir = tempfile.TemporaryDirectory(prefix=gff_name + ".", dir=self.tempdir)
            with tempfile.NamedTemporaryFile(
                                             prefix=gff_name + ".", 
                                             suffix=".agat", 
                                             dir=reading_temp_dir.name,
                                             delete=True
                                            ) as tmp:

                if Path(tmp.name).is_file():
                    os.remove(tmp.name)

                filename = Path(tmp.name).name
                command = f"agat_convert_sp_gxf2gxf.pl -g {gff} -o {filename}"

                # move inside tempdir
                cur_dir = os.getcwd()
                os.chdir(reading_temp_dir.name)

                subprocess.run(
                           command, 
                           check=True, 
                           shell=True, 
                           stderr=subprocess.DEVNULL, 
                           stdout=subprocess.DEVNULL
                        )
                
                gff_df = self.read_gff(filename)

            os.chdir(cur_dir)
            reading_temp_dir.cleanup()

        else:
            gff_df = self.read_gff(gff)


        if gff_df.shape[0] == 0:

            # read it again
            gff_df = self.read_gff(gff)
            compartments = gff_df['compartment'].unique()
            if len(compartments) > 1:
                raise ValueError("Invalid gff compartments after AGAT filtering. Must contain at least 'region'. Accession: {gff}.")

            if gff_df.shape[0] == 0:
                raise ValueError("Empty gff dataframe for accession {gff}.")


        gff_df.loc[:, "start"] = gff_df["start"] - 1
        gff_df.loc[:, "biotype"] = gff_df["attributes"].apply(GFFCleaner.parse_biotype)

        if reduce_biotype:
            gff_df.loc[:, "biotype"] = gff_df["biotype"].apply(GFFCleaner.reduce_biotype)

        gff_df = self.merge_compartments(gff_df)
        pybedtools.helpers.cleanup(verbose=False, remove_all=False)

        return gff_df

    def merge_compartments(self, gff_df: pd.DataFrame) -> pd.DataFrame:
        merged_gff = []
        merge_columns = ["seqID", "start", "end", "compartment", "source"]
        for compartment, biotype in self.all_compartments:

            if compartment == "gene":
                temp_comp = gff_df[(gff_df["compartment"] == compartment) | (gff_df["compartment"] == "pseudogene")]
                temp_comp.loc[:, "compartment"] = temp_comp["compartment"].replace("pseudogene", compartment)
            else:
                temp_comp = gff_df[gff_df["compartment"] == compartment]

            if biotype is not None:
                temp_comp = temp_comp[temp_comp["biotype"] == biotype]

            if temp_comp.shape[0] == 0:
                continue

            # temp_comp.reset_index(drop=True, inplace=True)
            merged_tool = BedTool.from_dataframe(temp_comp[merge_columns])\
                            .sort()\
                            .merge(
                                    c=[4, 5], 
                                    o=["count", "distinct"]
                            )

            merged_gff_comp = pd.read_table(
                                merged_tool.fn,
                                header=None,
                                names=["seqID", "start", "end", "overlapCount", "source"]
                                )
            merged_gff_comp.loc[:, "compartment"] = compartment 
            merged_gff_comp.loc[:, "biotype"] = biotype if biotype else "."
            merged_gff.append(merged_gff_comp)
    
        if len(merged_gff) > 0:
            merged_gff = pd.concat(merged_gff, axis=0)
            merged_gff.loc[:, "compartment"] = pd.Categorical(merged_gff["compartment"], categories=["region", "gene", "exon", "CDS"], ordered=True)
            
            merged_gff = merged_gff.sort_values(by=['seqID', 'start', 'compartment'], ascending=True)\
                                                      .reset_index(drop=True)
            merged_gff.loc[:, "compartment"] = merged_gff["compartment"].astype(str)

        else: 
            merged_gff = pd.DataFrame([], 
                                      columns=[
                                            "seqID", 
                                            "start", 
                                            "end", 
                                            "overlapCount",
                                            "source",
                                            "compartment", 
                                            "biotype"
                                            ]
                                      )

        return merged_gff


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=str)
    parser.add_argument("--add_exons", type=int, default=1)
    parser.add_argument("--bedtools_path", type=str, default="/storage/group/izg5139/default/nicole/miniconda3/bin")
    parser.add_argument("--tempdir", type=str, default="gff_clean_tmp")


    args = parser.parse_args()
    gff = args.gff
    add_exons = args.add_exons
    bedtools_path = args.bedtools_path
    tempdir_path = args.tempdir

    # bedtools_path = "/home/dollzeta/frogtools/bedtools2/bin"

    cleaner = GFFCleaner(bedtools_path=bedtools_path, tempdir=tempdir_path)
    gff_df = cleaner.read(gff, add_exons=add_exons)

    breakpoint()
