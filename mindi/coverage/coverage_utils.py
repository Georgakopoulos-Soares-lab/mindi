import os
import pandas as pd
from mindi.coverage.gff_clean import GFFCleaner
from pathlib import Path
from typing import Optional
from pybedtools import BedTool

COVERAGE_FIELDS = [
                   "seqID", 
                   "start", 
                   "end", 
                   "compartment", 
                   "biotype",
                   "overlapCount", 
                   "totalHits", 
                   "bpCovering", 
                   "compartmentLength", 
                   "coverage"
                   ]
extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
generate_columns = lambda partition_col: [
                          partition_col,
                          "compartment",
                          "biotype",
                          "compartmentLength",
                          "totalCompartments",
                          "bpCovering",
                          "atLeastOne",
                          "totalHits",
                          "totalCoverage",
                          "overlapping",
                          "averageCoverage",
                          "medianCoverage",
                          "minCoverage",
                          "maxCoverage"
                        ]
def extract_coverage(extraction_file: str | os.PathLike[str], 
                     gff_file: str | os.PathLike[str],
                     partition_col_collection: list[str],
                     partition_col: Optional[str] = None,
                     merge_compartments: bool = True,
                     is_merged: bool = False,
                     partition_on_biotype: bool = False) -> pd.DataFrame:
    global generate_columns, extract_id
    coverage_columns = generate_columns(partition_col)
    usecols = ["seqID", "start", "end"]
    if partition_col is not None:
        usecols.append(partition_col)
    extract_df = pd.read_table(extraction_file, usecols=usecols)
    gff_cleaner = GFFCleaner()
    gff_df_merged = gff_cleaner.read(gff_file, 
                                     add_exons=False, 
                                     merge_compartments=merge_compartments,
                                     is_merged=is_merged,
                                     partition_on_biotype=partition_on_biotype
                                     )
    coverage_table = []
    if gff_df_merged.shape[0] == 0:
        print(f"Accession {accession} was found without any relevant genomic compartments. Skipping...")
        coverage_table = pd.DataFrame([], columns=coverage_columns)
        return coverage_table
    accession_id = extract_id(extraction_file)
    # partition_col_collection = [0, 1, 2, 3, 4, 5, 6, 7, 8, all]
    for split_value in partition_col_collection:
        if split_value == "all":
            temp = extract_df
        else:
            temp = extract_df[extract_df[partition_col] == int(split_value)]
        extract_df_temp = BedTool.from_dataframe(temp)
        compartment_bed = (
                            BedTool.from_dataframe(gff_df_merged[["seqID", 
                                                                "start", 
                                                                "end", 
                                                                "compartment", 
                                                                "biotype", 
                                                                "overlapCount"
                                                            ]]
                                                   )
                                    .sort()
                        )
        coverage_df = pd.read_table(
                                compartment_bed.coverage(extract_df_temp).fn,
                                header=None,
                                names=COVERAGE_FIELDS)
        coverage_df.loc[:, "hadHit"] = (coverage_df["totalHits"] > 0).astype(int)
        coverage_df.loc[:, "coverage"] = 1e6 * coverage_df["coverage"]
        coverage_df = (
                    coverage_df.groupby(["compartment", "biotype"], as_index=False)
                                        .agg(
                                            totalHits=("totalHits", "sum"),
                                            atLeastOne=("hadHit", "sum"),
                                            bpCovering=("bpCovering", "sum"),
                                            compartmentLength=("compartmentLength", "sum"),
                                            totalCompartments=("overlapCount", "sum"),
                                            averageCoverage=("coverage", "mean"),
                                            medianCoverage=("coverage", "median"),
                                            minCoverage=("coverage", "min"),
                                            maxCoverage=("coverage", "max"))
                        )
        coverage_df.loc[:, partition_col] = split_value
        # coverage_df.loc[:, "compartment"] = coverage_df["compartment"].replace("region", "Genome")
        coverage_df.loc[:, "totalCoverage"] = 1e6 * coverage_df["bpCovering"].div(coverage_df["compartmentLength"])
        coverage_df.loc[:, "overlapping"] = 1e2 * coverage_df["atLeastOne"].div(coverage_df["totalCompartments"])
        coverage_df.loc[:, "#assembly_accession"] = accession_id 
        round_cols = ["averageCoverage", 
                      "medianCoverage", 
                      "minCoverage", 
                      "maxCoverage", 
                      "totalCoverage", 
                      "overlapping"]
        coverage_df[round_cols] = coverage_df[round_cols].round(3)
        coverage_df.set_index("#assembly_accession", inplace=True)
        coverage_table.append(coverage_df)
    if len(coverage_table) > 0:
        coverage_table = pd.concat(coverage_table, axis=0)
        coverage_table = coverage_table[coverage_columns]
    else:
        coverage_table = pd.DataFrame([], columns=coverage_columns)
    return coverage_table

if __name__ == "__main__":
    import random
    extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
    parent = Path(__file__).parent.parent
    extracted_file = random.choice([file for file in parent.joinpath("extractions").glob("*.tsv") if "IR" in file.name])
    print(extracted_file)
    extracted_df = pd.read_table(extracted_file)
    gff_files = {file.name.split('.gff')[0]: file for file in parent.joinpath("accessions").glob("*.gff")}
    chosen_file_id = extract_id(extracted_file)
    associated_gff = gff_files[chosen_file_id]
    coverage_df = extract_coverage(extracted_file, 
                                   associated_gff, 
                                   partition_col=None, # "spacerLength",
                                   partition_col_collection=["all"], # list(map(str, range(9))) + ["all"]
                                   )
