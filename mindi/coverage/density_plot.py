import pandas as pd
from pybedtools import BedTool
from mindi.coverage.gff_clean import GFFCleaner
from mindi.coverage.pwm_density import PWMExtractor, strand_evaluators
from mindi.coverage.windows_maker import WindowMaker
from dataclasses import dataclass 
from typing import Iterator, Optional, Callable
from pathlib import Path

@dataclass(frozen=True, slots=True, kw_only=True)
class Density:
    loci: str
    density: list[int]
    category: Optional[str] = None

FIELDS = ["seqID", 
          "start", 
          "end", 
          "strand", 
          "chromosome", 
          "motif_start", 
          "motif_end", 
          "sequence", 
          "overlap"]

def relative_density(center_extraction: str, extraction: str):
    maker = WindowMaker(base=1)
    extractions_origin_bed = BedTool.from_dataframe(
                            pd.read_table(extraction,
                                          usecols=["seqID", "start", "end", "sequence"]
                                        )
                            )
    extractions_bed = BedTool.from_dataframe(
                            pd.read_table(extraction,
                                          usecols=["seqID", "start", "end", "sequence"]
                                        )
                            )
    
def extract_density(extraction: str,
                    gff_file: str,
                    window_size: int,
                    density_type: str = "density",
                    determine_strand: Optional[Callable[[str], str]] = None,
                    genome: Optional[str] = None,
                    ) -> Iterator[Density]:
    global FIELDS
    maker = WindowMaker(base=1)
    cleaner = GFFCleaner(valid_compartments=["gene"])
    extractions_bed = BedTool.from_dataframe(
                            pd.read_table(extraction,
                                          usecols=["seqID", "start", "end", "sequence"]
                                        )
                            )
    gff_df = cleaner.read_gff(gff_file)[["seqID", "start", "end", "strand"]]
    gff_tss = maker.make_windows(gff_df, loci="start", genome=genome)
    gff_tes = maker.make_windows(gff_df, loci="end", genome=genome)
    gff_tss_bed = BedTool.from_dataframe(gff_tss)
    gff_tes_bed = BedTool.from_dataframe(gff_tes)
    intersect_df_tss = pd.read_table(
                        gff_tss_bed.intersect(extractions_bed, wo=True).fn,
                        header=None,
                        names=FIELDS
                        )
    intersect_df_tes = pd.read_table(
                            gff_tes_bed.intersect(extractions_bed, wo=True).fn,
                            header=None,
                            names=FIELDS
                            )
    pwm = PWMExtractor()
    if density_type == "template":
        density_evaluator = pwm.extract_template_density
    elif density_type == "pwm":
        density_evaluator = pwm.extract_PWM
    else:
        density_evaluator = pwm.extract_density

    for loci, intersect_df in zip(["tss", "tes"], [intersect_df_tss, intersect_df_tes]):
        if density_type == "template":
            density = density_evaluator(intersect_df, 
                                        window_size=window_size, 
                                        determine_strand=determine_strand)
        else:
            density = density_evaluator(intersect_df, 
                                        window_size=window_size)
        
        if isinstance(density, dict):
            for key, value in density.items():
                yield Density(density=value, loci=loci, category=key)
        else:
            yield Density(density=density, loci=loci)

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        window_size = int(sys.argv[1])
    else:
        window_size = 500
    density_type = "density"
    determine_strand = strand_evaluators["HDNA"]
    extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
    parent = Path(__file__).parent.parent
    extracted_file = [file for file in parent.joinpath("extractions").glob("*.tsv") if "IR" in file.name and "GCF_004355385.1" in file.name][0]
    print(extracted_file)
    gff_files = {file.name.split('.gff')[0]: file for file in parent.joinpath("accessions").glob("*.gff")}
    chosen_file_id = extract_id(extracted_file)
    associated_gff = gff_files[chosen_file_id]
    densities = dict()
    # if the name is the same it will overwrite the class above watch out?
    for eDensity in extract_density(
            extraction=extracted_file,
            gff_file=associated_gff,
            window_size=window_size,
            density_type=density_type,
                              determine_strand=determine_strand,
                              ):
        densities.update({(eDensity.loci, eDensity.category): eDensity.density}) 
