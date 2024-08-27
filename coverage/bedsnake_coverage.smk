# Coverage PIPELINE

__author__ = "Nikol Chantzi"
__email__ = "nmc6088@psu.edu"
__version__ = "1.0.1"

### > IMPORTS BEGIN

from pathlib import Path
import subprocess
from termcolor import colored
import threading
import os
import tempfile
import json
import pybedtools
from utils import ProgressTracker
from pybedtools import BedTool
import numpy as np
from gff_clean import GFFCleaner
import pandas as pd
from scheduling import MiniBucketScheduler

### > IMPORTS END

out = Path(config['out']).resolve()
out.mkdir(exist_ok=True)
TOTAL_BUCKETS = int(config['buckets'])
mode = config['mode']

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

rule schedule:
    output:
        '%s/schedule_coverage_%s.json' % (out, TOTAL_BUCKETS)
    params:
        files=config['files'],
        gff_parent=Path(config['gff_parent_coverage']).resolve(),
    run:
        assemblies = []
        with open(params.files, mode='r', encoding='UTF-8') as f:
            for line in f:
                line = line.strip()

                if line.count("\t") > 0:
                    line = line.split("\t")[0]

                # is associated to GFF
                gff_corresponding = params.gff_parent.joinpath(Path(line).name.replace("fna.gz", "agat.merged.gff"))

                if gff_corresponding.is_file():
                    assemblies.append(str(gff_corresponding))

        if len(assemblies) == 0:
            color = "red"
        else:
            color = "green"

        print(colored(f"Total assemblies detected: {len(assemblies)}.", color))
        if len(assemblies) == 0:
            raise ValueError(f'No assemblies were detected from the path {params.files}.')

        mini_bucket_scheduler = MiniBucketScheduler()
        scheduled_files = mini_bucket_scheduler.schedule(assemblies, total_buckets=TOTAL_BUCKETS)
        mini_bucket_scheduler.saveas(scheduled_files, output[0])


GFF_FIELDS = ["seqID", "source", "compartment", "start", "end", "score", "strand", "phase", "attributes"]
COVERAGE_FIELDS = ["seqID", "start", "end", "compartment", "biotype", "overlapCount", "totalHits", "bpCovering", "compartmentLength", "coverage"]

rule mergeGFF:
    input:
        '%s/schedule_coverage_%s.json' % (out, TOTAL_BUCKETS)
    output:
        '%s/flow/bucket_{bucket}.%s.merged.completed' % (out, mode)
    params:
        tempdir=Path(config['tempdir']).resolve(),
        bedtools_path=Path(config['bedtools_path']).resolve(),
        destination=Path(config['destination_merged']).resolve()
    run:
        Path(f"{out}/flow").mkdir(exist_ok=True)
        all_compartments = [
                            ("region", None),
                            ("gene", None),
                            ("five_prime_UTR", None),
                            ("three_prime_UTR", None),
                            ("gene", "protein_coding"),
                            ("gene", "non_coding"),
                            ("exon", None),
                            ("CDS", None),
            ]
        gff_cleaner = GFFCleaner(
                            tempdir=params.tempdir,
                            bedtools_path=params.bedtools_path,
                            all_compartments=all_compartments,
                        )

        buckets = load_bucket(wildcards.bucket)
        def extract_name(accession):
            accession = Path(accession).name
            return accession.split(".agat.gff")[0]

        for accession in buckets:
            accession_name = extract_name(accession)
            destination = params.destination.joinpath(accession_name + ".agat.merged.gff")

            if destination.is_file():
                continue

            merged_gff = gff_cleaner.read(accession, add_exons=False)
            merged_gff.to_csv(destination, mode="w", index=False, sep="\t", header=True)


rule extractCoverage:
    input:
        '%s/schedule_coverage_%s.json' % (out, TOTAL_BUCKETS)
    output:
        '%s/%s/coverage/coverage_compartments_bucket_{bucket}.%s.coverage' % (out, mode, mode)
    params:
        out=Path(config['out']).resolve(),
        gff_parent=Path(config['gff_parent_coverage']).resolve(),
        extraction_parent=Path(config['extraction_parent']).resolve(),
        tempdir=Path(config['tempdir']).resolve(),
        bedtools_path=config['bedtools_path'],
        split_category=config['split_category'],
        split_collection=config['split_collection'],
        mode=config["mode"],
    run:
        extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])

        path_to_glob = params.extraction_parent.joinpath(f'{params.mode}_extracted_accessions')
        extractions = {extract_id(file): file for file in path_to_glob.glob("*.csv")}

        accessions = load_bucket(wildcards.bucket)
        util_cols = ["seqID", "start", "end"]

        print(colored(f"Splitting coverage process on column {params.split_category}.", "blue"))
        split_category_collection = list(map(str, params.split_collection)) + ["all"]
        coverage_table = []

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

        for accession in accessions:

            tracker.counter += 1
            print(colored(f"Processing accession '{accession}'.", "green"))

            accession = Path(accession)
            accession_id = extract_id(accession)
            gff_file = accession

            # gff_file = params.gff_parent.joinpath(accession.name.split(".gff")[0] + ".agat.merged.gff")
            if accession_id not in extractions:
                print(f"Accession {accession_id} has not been extracted. Skipping...")
                continue

            extraction_file = extractions[accession_id]

            # if not gff_file.is_file():
            #    print(f"AGAT {gff_file}.")
            #    raise FileNotFoundError(f'Could not locate gff file for accession {accession}')

            # extraction file
            extract_df = pd.read_table(
                                       extraction_file,
                                       usecols=["seqID", "start", "end", params.split_category]
                                    )

            # accession_name = extract_name(accession)
            gff_df_merged = pd.read_table(
                                    gff_file,
                                    header=None,
                                    names=["seqID", "start", "end", "overlapCount", "source", "compartment", "biotype"],
                                    dtype={
                                            "start": int,
                                            "end": int
                                        }
                                    )
            if gff_df_merged.shape[0] == 0:
                print(f"Accession {accession} was found without any relevant genomic compartments. Skipping...")
                continue

            # split_category_collection = [0, 1, 2, 3, 4, 5, 6, 7, 8, all]

            for split_value in split_category_collection:
                if split_value == "all":
                    temp = extract_df
                else:
                    split_value = int(split_value)
                    temp = extract_df[extract_df[params.split_category] == split_value]

                extract_df_temp = BedTool.from_dataframe(temp)
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
                                    )

                coverage_df.loc[:, params.split_category] = split_value

                coverage_df.loc[:, "compartment"] = coverage_df["compartment"].replace("region", "Genome")
                coverage_df.loc[:, "totalCoverage"] = 1e6 * coverage_df["bpCovering"].div(coverage_df["compartmentLength"])
                coverage_df.loc[:, "overlapping"] = 1e2 * coverage_df["atLeastOne"].div(coverage_df["totalCompartments"])
                coverage_df.loc[:, "#assembly_accession"] = extract_id(accession)
                round_cols = ["averageCoverage", "medianCoverage", "minCoverage", "maxCoverage", "totalCoverage", "overlapping"]
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

