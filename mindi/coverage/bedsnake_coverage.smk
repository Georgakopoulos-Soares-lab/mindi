# Coverage PIPELINE

__author__ = "Nikol Chantzi"
__email__ = "nmc6088@psu.edu"
__version__ = "1.0.1"

### > IMPORTS BEGIN

from pathlib import Path
import subprocess
from termcolor import colored
import csv
import threading
import os
import tempfile
import json
import pybedtools
import shutil
from utils import ProgressTracker
from pybedtools import BedTool
import numpy as np
import pandas as pd
from mindi.coverage.gff_clean import GFFCleaner
# from mindi.coverage.density import COVERAGE_FIELDS, GFF_FIELDS
from mindi.scheduling import MiniBucketScheduler

### > IMPORTS END
GFF_FIELDS = ["seqID", "source", "compartment", "start", "end", "score", "strand", "phase", "attributes"]
COVERAGE_FIELDS = ["seqID", "start", "end", "compartment", "biotype", "overlapCount", "totalHits", "bpCovering", "compartmentLength", "coverage"]

out = Path(config['out']).resolve()
out.mkdir(exist_ok=True)
TOTAL_BUCKETS = int(config['buckets'])
mode = config['mode']
DESIGN = config['DESIGN']

tempdir = Path(config['tempdir']).resolve()
tempdir.mkdir(exist_ok=True)
pybedtools.set_tempdir(config['tempdir'])

def load_bucket(bucket) -> list[str]:
    global TOTAL_BUCKETS
    global out
    with open(f"{out}/schedule_coverage_{TOTAL_BUCKETS}.json", mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket)]

extract_id = lambda assembly: '_'.join(Path(assembly).name.split('_')[:2])

rule all:
    input:
        '%s/%s/coverage/coverage_compartments.%s.parquet' % (out, mode, mode)

def get_file_ids(design: str) -> dict[str, str]:
    file_ids = {}
    with open(design, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            file_ids.update({row['accession_id']: row['extraction']})
    return file_ids

rule schedule:
    input:
      DESIGN
    output:
        '%s/schedule_coverage_%s.json' % (out, TOTAL_BUCKETS)
    run:
        assemblies = []
        accession_ids = []
        with open(input[0], mode='r', encoding='UTF-8') as f:
            reader = csv.DictReader(f, delimiter=",")
            for row in reader:
              assemblies.append(row['gff'])
        print(colored(f"Total assemblies detected: {len(assemblies)}.", "green"))
        if len(assemblies) == 0:
            raise ValueError(f'No assemblies were detected from the path {params.files}.')
        scheduler = MiniBucketScheduler()
        scheduled_files = scheduler.schedule(assemblies, total_buckets=TOTAL_BUCKETS)
        scheduler.saveas(scheduled_files, output[0])

rule extractCoverage:
    input:
        DESIGN, '%s/schedule_coverage_%s.json' % (out, TOTAL_BUCKETS)
    output:
        '%s/%s/coverage/coverage_compartments_bucket_{bucket}.%s.coverage' % (out, mode, mode)
    params:
        out=Path(config['out']).resolve(),
        tempdir=Path(config['tempdir']).resolve(),
        split_category=config['split_category'],
        split_collection=config['split_collection'],
        mode=config["mode"],
        # bedtools_path=config['bedtools_path'],
    run:
        extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])
        extract_gff = lambda gff: extract_id(Path(gff).name.split('.gff')[0])
        accessions = load_bucket(wildcards.bucket)
        file_ids = get_file_ids(input[0])

        all_compartments = [
                                 ("region", None),
                                 ("CDS", None),
                                 ("exon", None),
                                 ("gene", None),
                                 ("five_prime_UTR", None),
                                 ("three_prime_UTR", None),
                                 ("gene", "protein_coding"),
                                 ("gene", "non_coding")
                            ]

        print(colored(f"Splitting coverage process on column {params.split_category}.", "blue"))
        partition_col = ["all"]
        if len(params.split_collection) > 0:
          partition_col = list(map(str, params.split_collection))
        
        # >> LOGGING START
        total_accessions = len(accessions)
        progress_tracking_log = params.out.joinpath("biologs", f"biolog_tracker_{mode}_{wildcards.bucket}.coverage.log")
        tracker = ProgressTracker(
                                 total_accessions=total_accessions,
                                 filename=progress_tracking_log,
                                 bucket_id=wildcards.bucket,
                                 sleeping_time=300
                                 )
        logging_thread = threading.Thread(target=tracker.track_progress, daemon=True)
        logging_thread.start()
        # << LOGGING END

        gff_reader = GFFCleaner(all_compartments=all_compartments)
        coverage_table = []
        for gff_file in accessions:
            tracker.counter += 1
            print(colored(f"Processing accession '{gff_file}'.", "green"))

            accession_id = extract_gff(gff_file)
            extraction_file = file_ids[accession_id]

            extract_df = pd.read_table(extraction_file) #, usecols=["seqID", "start", "end", params.split_category])
            # accession_name = extract_name(accession)
            gff_df_merged = gff_reader.read(gff_file, 
                                            add_exons=False,
                                            change_names=False,
                                            merge_compartments=True,
                                            partition_on_biotype=True,
                                            is_merged=False)
            if gff_df_merged.shape[0] == 0:
                print(colored(f"Accession {accession} was found without any relevant genomic compartments. Skipping...", "red"))
                continue

            # partition_col = [0, 1, 2, 3, 4, 5, 6, 7, 8, all]
            for split_value in partition_col:
                if split_value == "all":
                    temp_df = extract_df
                else:
                    if split_value.isdigit():
                        split_value = int(split_value)
                    temp_df = extract_df[extract_df[params.split_category] == split_value]
                extract_df_temp = BedTool.from_dataframe(temp_df).sort()
                compartment_df = BedTool.from_dataframe(gff_df_merged[[
                                                                       "seqID",
                                                                       "start",
                                                                       "end",
                                                                       "compartment",
                                                                       "biotype",
                                                                       "overlapCount"
                                                                       ]]
                                                        )\
                                                .sort()
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
                                            stdCoverage=("coverage", "std"),
                                    )
                coverage_df.loc[:, params.split_category] = split_value
                coverage_df.loc[:, "compartment"] = coverage_df["compartment"].replace("region", "Genome")
                coverage_df.loc[:, "totalCoverage"] = 1e6 * coverage_df["bpCovering"].div(coverage_df["compartmentLength"])
                coverage_df.loc[:, "overlapping"] = 1e2 * coverage_df["atLeastOne"].div(coverage_df["totalCompartments"])
                coverage_df.loc[:, "#assembly_accession"] = extract_id(gff_file)
                round_cols = ["averageCoverage", "medianCoverage", "minCoverage", "maxCoverage", "totalCoverage", "overlapping", "stdCoverage"]
                coverage_df[round_cols] = coverage_df[round_cols].round(3)
                coverage_table.append(coverage_df)

        coverage_columns = [
                          "#assembly_accession",
                          params.split_category,
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
        if len(coverage_table) > 0:
            coverage_table = pd.concat(coverage_table, axis=0)
            coverage_table = coverage_table[coverage_columns]
        else:
            coverage_table = pd.DataFrame([], columns=coverage_columns)
        pybedtools.helpers.cleanup(remove_all=False)
        coverage_table.set_index("#assembly_accession", inplace=True)
        coverage_table.to_csv(output[0], sep="\t", index=True, mode="w")


rule mergeCoverage:
    input:
        expand('%s/%s/coverage/coverage_compartments_bucket_{bucket}.%s.coverage' % (out, mode, mode),
               bucket=range(TOTAL_BUCKETS))
    output:
        '%s/%s/coverage/coverage_compartments.%s.parquet' % (out, mode, mode)
    run:
        coverage_table = []
        for bucket in range(TOTAL_BUCKETS):
            coverage_file = f"{out}/{mode}/coverage/coverage_compartments_bucket_{bucket}.{mode}.coverage"
            coverage_df = pd.read_table(coverage_file)
            coverage_table.append(coverage_df)
        coverage_table = pd.concat(coverage_table, axis=0)
        coverage_table.to_parquet(output[0], engine="fastparquet")
