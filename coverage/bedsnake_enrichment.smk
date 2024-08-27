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
from pybedtools import BedTool
from utils import ProgressTracker
import numpy as np
import pandas as pd
from pwm_density import PWMExtractor
from collections import defaultdict
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
    with open(f"{out}/schedule_enrichment_{TOTAL_BUCKETS}.json", mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket)]

extract_id = lambda assembly: '_'.join(Path(assembly).name.split('_')[:2])

rule all:
    input:
        '%s/%s/enrichment/enrichment_compartments.TSS.%s.parquet' % (out, mode, mode),
        '%s/%s/enrichment/enrichment_compartments.TES.%s.parquet' % (out, mode, mode)

rule schedule:
    output:
        '%s/schedule_enrichment_%s.json' % (out, TOTAL_BUCKETS)
    params:
        files=config['files'],
        gff_parent=Path(config['gff_parent_enrichment']).resolve(),
    run:
        assemblies = []
        with open(params.files, mode='r', encoding='UTF-8') as f:
            for line in f:
                line = line.strip()

                if line.count("\t") > 0:
                    line = line.split("\t")[0]

                # is associated to GFF
                gff_corresponding = params.gff_parent.joinpath(Path(line).name.replace("fna", "gff"))

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


GFF_FIELDS = [
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

INTERSECT_FIELDS = ["seqID",
                    "start",
                    "end",
                    "strand",
                    "biotype",
                    "chromosome",
                    "motif_start",
                    "motif_end",
                    "sequenceOfArm",
                    "spacerLength",
                    "overlap"]


def parse_biotype(attributes: str) -> str:
    if "gene_biotype" in attributes:
        gene_biotype = attributes.split("gene_biotype=")[1].split(";")[0]

        if gene_biotype == "protein_coding":
            return gene_biotype

        return "non_coding"

    return "."

def extract_parent_id(attributes: str) -> str:
    return ""

rule extractEnrichment:
    input:
        '%s/schedule_enrichment_%s.json' % (out, TOTAL_BUCKETS)
    output:
        '%s/%s/enrichment/enrichment_compartments_bucket_{bucket}.%s.enrichment' % (out, mode, mode)
    params:
        out=Path(config['out']).resolve(),
        extraction_parent=Path(config['extraction_parent']).resolve(),
        tempdir=Path(config['tempdir']).resolve(),
        bedtools_path=config['bedtools_path'],
        compartment=config['compartment'],
        split_category=config['split_category'],
        split_collection=config['split_collection'],
        window_size=int(config["window_size"]),
        mode=config["mode"],
    run:
        print(f"Initializing bucket {wildcards.bucket} enrichment extraction process for compartment {params.compartment} and window length {params.window_size}.")
        extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])
        pybedtools.set_bedtools_path(path=params.bedtools_path)
        pybedtools.helpers.set_bedtools_path(path=params.bedtools_path)

        path_to_glob = params.extraction_parent.joinpath(f'{params.mode}_extracted_accessions')
        extractions = {extract_id(file): file for file in path_to_glob.glob("*.csv")}

        accessions = load_bucket(wildcards.bucket)
        util_cols = ["seqID", "start", "end"]

        print(colored(f"Splitting coverage process on column {params.split_category}.", "blue"))
        split_category_collection = list(map(int, params.split_collection))
        enrichment_table = []

        total_accessions = len(accessions)
        progress_tracking_log = params.out.joinpath("biologs", f"biolog_tracker_{mode}_{wildcards.bucket}.enrichment.log")
        tracker = ProgressTracker(
                                 total_accessions=total_accessions,
                                 filename=progress_tracking_log,
                                 bucket_id=wildcards.bucket,
                                 sleeping_time=300
                                 )

        logging_thread = threading.Thread(target=tracker.track_progress, daemon=True)
        logging_thread.start()

        transcription_sites = ["TSS", "TES"]
        pwm = PWMExtractor(mode=params.mode)

        if mode == "STR":
            sequence_col = "sequence"
        elif mode == "DR" or mode == "MR" or mode == "IR":
            sequence_col = "sequenceOfArm"
        else:
            raise ValueError(f'Unknown mode {mode}.')

        for gff_file in accessions:
            tracker.counter += 1
            print(colored(f"Processing accession '{gff_file}'.", "green"))

            accession_id = extract_id(gff_file)
            if accession_id not in extractions:
                print(f"Accession {accession_id} has not been extracted. Skipping...")
                continue

            extraction_file = extractions[accession_id]

            # extraction file
            extract_df = pd.read_table(
                                       extraction_file,
                                       usecols=["seqID", "start", "end", sequence_col, params.split_category, "sequence"]
                                    )
            if params.mode == "MR":

                print("Filtering H-DNA!")
                # H-DNA filtering
                extract_df.loc[:, "sequenceLength"]  = extract_df["sequence"].apply(len)
                extract_df.loc[:, "at_content"] = extract_df["sequence"].str.count("a|t")
                extract_df.loc[:, "at_content"] = extract_df["at_content"].div(extract_df["sequenceLength"])

                extract_df.loc[:, "pyrine"] = extract_df["sequence"].str.count("g|a").div(extract_df["sequenceLength"])
                extract_df.loc[:, "pyrimidine"] = extract_df["sequence"].str.count("c|t").div(extract_df["sequenceLength"])

                extract_df = extract_df[((extract_df["pyrine"] >= 0.9) | (extract_df["pyrimidine"] >= 0.9)) & (extract_df["at_content"] <= 0.8)]

                if extract_df.shape[0] == 0:
                    print(f"No H-DNA for accession {gff_file}! Skipping...")
                    continue

                extract_df = extract_df.reset_index(drop=True)\
                                        .drop(columns=[
                                                "pyrine",
                                                "pyrimidine",
                                                "at_content",
                                                "sequenceLength",
                                                "sequence",
                                                ]
                                        )


            # accession_name = extract_name(accession)
            gff_df = pd.read_table(
                                    gff_file,
                                    header=None,
                                    comment="#",
                                    names=GFF_FIELDS,
                                    dtype={
                                            "start": int,
                                            "end": int
                                        }
                            )

            unique_compartments = set(gff_df.compartment)
            if params.compartment not in unique_compartments:
                print(colored(f"Invalid compartment detected {params.compartment} for accession '{gff_file}'.", "red"))
                continue

            gff_df = gff_df[gff_df["compartment"] == params.compartment]\
                                .reset_index(drop=True)\
                                .drop(columns=["compartment"])
            gff_df.loc[:, "biotype"] = gff_df["attributes"].apply(parse_biotype)
            gff_df.loc[:, "start"] = gff_df["start"] - 1
            gff_df.loc[:, "end"] = gff_df["end"] - 1

            print(colored(f"Total {gff_df.shape[0]} {params.compartment}(s) detected for accession '{gff_file}'.", "green"))

            # filter first and last exon
            if params.compartment == "exon":
                gff_df.loc[:, "parentID"] = gff_df["attributes"].apply(extract_parent_id)


            gene_biotypes = set(gff_df["biotype"])
            def _fetch_origin(row, site):

                if site != 'TSS' and site != 'TES':
                    raise ValueError(f'Unknown transcription site {site}.')

                # At the GENIC START return:
                # the start of GFF when strand is positive
                # the end of GFF when strand is negative
                if site == 'TSS':
                    if row['strand'] == '+':
                        return row['start']
                    return row['end']

                # At the GENIC END return:
                # the end of GFF when strand is positive
                # the start of GFF when strand is negative
                if row['strand'] == '+':
                    return row['end']
                return row['start']


            compartment_stats = []

            for site in transcription_sites:
                gff_df.loc[:, "origin"] = gff_df[["strand", "start", "end"]].apply(lambda row: _fetch_origin(row, site), axis=1)
                gff_df.loc[:, f"{site}_start"] = np.maximum(gff_df["origin"] - params.window_size, 0)
                gff_df.loc[:, f"{site}_end"] = gff_df["origin"] + params.window_size + 1

                extract_bed = BedTool.from_dataframe(extract_df)
                compartment_df = BedTool.from_dataframe(gff_df[[
                                                            "seqID",
                                                            f"{site}_start",
                                                            f"{site}_end",
                                                            "strand",
                                                            "biotype",
                                                    ]]
                                            )\
                                        .sort()

                intersect_df = pd.read_table(
                                        compartment_df.intersect(extract_bed, wo=True).fn,
                                        header=None,
                                        names=INTERSECT_FIELDS,
                                        dtype={
                                            "start": int,
                                            "end": int,
                                            "motif_start": int,
                                            "motif_end": int,
                                            params.split_category: int,
                                              }
                                        )

                # TS_aggregated_statistics = intersect_df.groupby(["seqID", "start", "end"], as_index=False)\
                #                                       .agg(
                #                                            overlap=("overlap", "sum"),
                #                                            totalGenes=("overlap", "count"),
                #                                        )
                # TS_aggregated_statistics.loc[:, "atLeastOne"] = (TS_aggregated_statistics["overlap"] > 0).astype(int)
                # comparment_stats.append(TS_aggregated_statistics["atLeastOne"].value_counts())
                # intersect_df = intersect_df[intersect_df["overlap"] > 0].reset_index(drop=True)

                for biotype in gene_biotypes:
                    for attribute in params.split_collection:

                        print(f"Bucket {wildcards.bucket}; Processing {biotype=};{params.split_category}={attribute}...")
                        intersect_temp = intersect_df[(intersect_df['biotype'] == biotype) & (intersect_df[params.split_category] == attribute)]
                        vector_counts = pwm.extract_PWM(intersect_temp, window_size=params.window_size)
                        vector_counts = pd.DataFrame(vector_counts).T\
                                .reset_index()\
                                .rename(columns={"index": "nucleotide"})\
                                .rename(columns={str(i): f"{int(i)-params.window_size}" for i in range(2*params.window_size+1)})
                        vector_counts.loc[:, "site"] = site
                        vector_counts.loc[:, params.split_category] = attribute
                        vector_counts.loc[:, "biotype"] = biotype
                        vector_counts.loc[:, "#assembly_accession"] = accession_id

                        prev_cols = set(vector_counts.columns)
                        # column rearrangement
                        vector_counts = vector_counts[["#assembly_accession",
                                                       "site",
                                                       "biotype",
                                                       params.split_category,
                                                       "nucleotide"] + list(range(-params.window_size, params.window_size+1))]

                        assert prev_cols == set(vector_counts.columns), "Invalid columns."
                        enrichment_table.append(vector_counts)

        if len(enrichment_table) > 0:
            enrichment_table = pd.concat(enrichment_table, axis=0)
        else:
            print(colored(f"No files detected for bucket {wildcards.bucket}.", "red"))
            enrichment_table = pd.DataFrame([], columns=["#assembly_accession",
                                                         "site",
                                                         "biotype",
                                                         params.split_category,
                                                         "nucleotide"] + list(range(-params.window_size, params.window_size+1))
        pybedtools.helpers.cleanup(remove_all=False)
        enrichment_table.set_index("#assembly_accession", inplace=True)
        enrichment_table.to_csv(output[0], sep="\t", index=True, mode="w")


rule reduceEnrichment:
    input:
        expand('%s/%s/enrichment/enrichment_compartments_bucket_{bucket}.%s.enrichment' % (out, mode, mode),
               bucket=range(TOTAL_BUCKETS))
    output:
        '%s/%s/enrichment/enrichment_compartments.TSS.%s.parquet' % (out, mode, mode),
        '%s/%s/enrichment/enrichment_compartments.TES.%s.parquet' % (out, mode, mode)
    run:
        enrichment_table = []
        for bucket in range(TOTAL_BUCKETS):
            coverage_file = f"{out}/{mode}/enrichment/enrichment_compartments_bucket_{bucket}.{mode}.enrichment"
            coverage_df = pd.read_table(coverage_file)
            enrichment_table.append(coverage_df)

        enrichment_table = pd.concat(enrichment_table, axis=0)
        enrichment_table.set_index("#assembly_accession", inplace=True)
        enrichment_table_TSS = enrichment_table[enrichment_table['site'] == 'TSS'].drop(columns=['site'])
        enrichment_table_TES = enrichment_table[enrichment_table['site'] == 'TES'].drop(columns=['site'])

        enrichment_table_TSS.to_parquet(output[0], engine="fastparquet")
        enrichment_table_TES.to_parquet(output[1], engine="fastparquet")


rule groupByEnrichment:
    input:
        '%s/%s/enrichment/enrichment_compartments.%s.parquet' % (out, mode, mode)
    output:
        '%s/%s/enrichment/enrichment_compartments.%s.grouped.parquet' % (out, mode, mode)
    params:
        window_size=int(config['window_size']),
        split_category=config['split_category'],
    run:
        enrichment_table = pd.read_parquet(input[0], engine="fastparquet")
        transcription_sites = ["start", "end"]

        for site in transcription_sites:
            grouped_enrichment = enrichment_table[enrichment_table["site"] == site]\
                                    .groupby(["biotype", "nucleotide", params.split_collection], as_index=False)\
                                    .agg({str(i): "sum" for i in range(-params.window_size, params.window_size+1)})
