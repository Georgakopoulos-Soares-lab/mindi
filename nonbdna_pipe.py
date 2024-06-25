from mindi import MindiHunter
from termcolor import colored
import json
import numpy as np
import pandas as pd
from utils import parse_fasta

buckets = config['buckets']
out = config['out']

def load_bucket(bucket):
    global buckets
    global out
    with open(f"{out}/schedule_{buckets}.json", mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket)]

extract_id = lambda assembly: '_'.join(Path(assembly).name.split('_')[:2])

def extract_name(accession):
    accession = Path(accession).name
    return accession.split(".fna")[0]

rule all:
    input:
        "%s/irp_completed/inverted_repeats_db.parquet.snappy" % out

rule schedule:
    output:
        "%s/schedule_%s.json" % (out, config["buckets"])
    params:
        buckets=config["buckets"],
        files=config["files"],
    run:
        print(colored(f"{len(params.files)} annotations detected.", "green"))
        assemblies = []
        with open(params.files, mode="r") as f:
            for line in f:
                line = line.strip()
                if line.count("\t") > 0:
                    line = line.split("\t")[0]

                assemblies.append(line)

        total_split = params.buckets
        splitted_batches = {batch_id: job.tolist() for batch_id, job in enumerate(np.array_split(assemblies, total_split))}
        
        Path(out).mkdir(exist_ok=True)
        with open(output[0], mode="w", encoding="UTF-8") as f:
            json.dump(splitted_batches, f, indent=4)


rule extractInvertedRepeats:
    input:
        '%s/schedule_%s.json' % (out, config["buckets"])
    output:
        touch("%s/irp_completed/bucket_{bucket}.completed" % out)
    params:
        out=Path(config['out']).resolve(),
        minrep=int(config['minrep']),
        max_spacer_length=int(config['max_spacer_length']),
        min_arm_length=int(config['min_arm_length']),
    run:
        Path(f"{params.out}/irp_completed").mkdir(exist_ok=True, parents=True)
        bucket = load_bucket(wildcards.bucket)

        tempdir = Path(params.out).joinpath("irp_temp")
        mindi = MindiHunter(tempdir=tempdir)

        destination_dir = params.out.joinpath("extracted_accessions")
        destination_dir.mkdir(exist_ok=True)

        for accession in bucket:
            accession_id = extract_name(accession)
            destination = destination_dir.joinpath(accession_id + ".IR.csv")

            _ = mindi.extract_inverted(
                          accession, 
                          minrep=params.minrep,
                          min_arm_length=params.min_arm_length,
                          max_spacer_length=params.max_spacer_length
                    )

            # .to_dataframe()

            irp_df = mindi.to_dataframe()
            total_chr = 0
            for seqID, sequence in parse_fasta(accession):

                temp = irp_df[irp_df['seqID'] == seqID]
                assert temp.shape[0] > 0
                total_chr = 1

                total = temp.shape[0]
                total_validated = 0

                for _, row in temp.iterrows():

                    start = row['start']
                    end = row['end']
                    derived_seq = row['sequence']
                    arm_seq = row['sequenceOfArm']
                    arm_len = len(arm_seq)
                    spacer_seq = row['sequenceOfSpacer']
                    if spacer_seq == ".":
                        spacer_seq = ""

                    spacer_true_len = row['spacerLength']
                    arm_true_len = row['armLength']
                    seq_true_len = row['sequenceLength']

                    spacer_len = len(spacer_seq)

                    section = sequence[start: end]
                    
                    assert section == derived_seq == arm_seq + spacer_seq + ''.join({'a': 't', 'g': 'c', 'c': 'g', 't': 'a'}[n] for n in arm_seq)[::-1]
                    assert arm_len * 2 + spacer_len == len(section) == len(derived_seq) == arm_true_len * 2 + spacer_true_len == seq_true_len
                    total_validated += 1

                assert total_validated == total

            assert total_chr == 1





            # irp_df.to_csv(
            #              destination,
            #              sep="\t",
            #              index=False
            #        )

            mindi.moveto(destination)

rule reduceInvertedRepeats:
    input:
        expand("%s/irp_completed/bucket_{bucket}.completed" % out, bucket=range(buckets))
    output:
        "%s/irp_completed/inverted_repeats_db.parquet.snappy" % out
    params:
        out=Path(config["out"]).resolve(),
        buckets=config["buckets"]
    run:
        files = [extracted_file for extracted_file in params.out.joinpath("extracted_accessions").glob("*.IR.csv")]
        df_all = []
        for file in files:
            df = pd.read_table(file)
            df_all.append(df)

        df_all = pd.concat(df_all, axis=0)
        df_all.to_parquet(output[0], engine="fastparquet")
