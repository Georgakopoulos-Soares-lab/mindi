import os
import pandas as pd
from gff_clean import GFFCleaner
import pybedtools
from pybedtools import BedTool


COVERAGE_FIELDS = ["seqID", "start", "end", "compartment", "counts", "totalHits", "bpCovering", "compartmentLength", "coverage"]

def extract_coverage(extraction_file: os.PathLike[str], 
                     gff_file: os.PathLike[str],
                     split_category: str,
                     split_category_collection: list[str]) -> pd.DataFrame:

    # extraction file
    extract_df = pd.read_table(
                            extraction_file, 
                            usecols=["seqID", "start", "end", split_category]
                    )
            
    # accession_name = extract_name(accession)
    gff_df_merged = gff_cleaner.read(gff_file, add_exons=False)
    coverage_table = []

    for split_value in split_category_collection:
        if split_value == "all":
            temp = extract_df
        else:
            split_value = int(split_value)
            temp = extract_df[extract_df[split_category] == split_value]

        extract_df_temp = BedTool.from_dataframe(temp)
        compartment_df = BedTool.from_dataframe(gff_df_merged[["seqID", 
                                                                "start", 
                                                                "end", 
                                                                "compartment", 
                                                                "biotype", 
                                                                "overlapCount"
                                                            ]]
                                                )

        coverage_df = pd.read_table(
                                compartment_df.coverage(extract_df_temp).fn,
                                header=None,
                                names=COVERAGE_FIELDS,
                            )
        coverage_df.loc[:, "hadHit"] = (coverage_df["totalHits"] > 0).astype(int)
        coverage_df.loc[:, "coverage"] = 1e6 * coverage_df["coverage"]

        coverage_df = coverage_df.groupby(["compartment", "biotype"], as_index=False)\
                                        .agg(
                                            totalHits=("totalHits", "sum"),
                                            atLeastOne=("hadHit", "sum"),
                                            bpCovering=("bpCovering", "sum"),
                                            compartmentLength=("compartmentLength", "sum"),
                                            totalCompartments=("overlapCount", "sum"),
                                            averageCoverage=("coverage", "mean"),
                                            medianCoverage=("coverage", "median"),
                                            minCoverage=("coverage", "min"),
                                            maxCoverage=("coverage", "max"),
                                    )

        coverage_df.loc[:, split_category] = split_value
        coverage_df.loc[:, "compartment"] = coverage_df["compartment"].replace("region", "Genome")
        coverage_df.loc[:, "totalCoverage"] = 1e6 * coverage_df["bpCovering"].div(coverage_df["compartmentLength"])
        coverage_df.loc[:, "overlapping"] = 1e2 * coverage_df["atLeastOne"].div(coverage_df["totalCompartments"])
        coverage_df.loc[:, "#assembly_accession"] = extract_id(accession)

        round_cols = ["averageCoverage", "medianCoverage", "minCoverage", "maxCoverage", "coverage", "overlapping"]
        coverage_df[round_cols] = coverage_df[round_cols].round(3)
        coverage_df.set_index("#assembly_accession", inplace=True)
        coverage_table.append(coverage_df)

    if len(coverage_table) > 0:
        coverage_table = pd.concat(coverage_table, axis=0)
    else:
        coverage_table = pd.DataFrame([])

    return coverage_table
