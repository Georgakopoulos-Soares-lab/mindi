from pathlib import Path
import subprocess
import os
import tempfile
import pybedtools
from gff_clean import GFFCleaner

pybedtools.set_tempdir(config['tempdir'])

out = Path(config['out']).resolve()
total_buckets = int(config['buckets'])
mode = config['mode']

def load_bucket(bucket) -> list[str]:
    global buckets
    global out
    with open(f"{out}/schedule_{buckets}.json", mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket)]

extract_id = lambda assembly: '_'.join(Path(assembly).name.split('_')[:2])

rule schedule:
    output:
        '%s/schedule_%s.json' % (out, total_buckets)
    params:
        files=config['files'],
        gff_parent=Path(config['gff_parent']).resolve(),
    run:
        assemblies = []
        with open(params.files, mode='r', encoding='UTF-8') as f:
            for line in f:
                line = line.strip()
                
                if line.count("\t") > 1:
                    line = line.split("\t")[0]
                
                # is associated to GFF
                if gff_parent.joinpath(line.replace("fna", "gff")).is_file():
                    assemblies.append(line)

        splitted_jobs = {bucket_id: job.tolist() for bucket_id, job in enumerate(np.array_split(assemblies, total_buckets), 0)}
        with open(output[0], mode='w', encoding='UTF-8') as f:
            json.dump(assemblies, f, indent=4)

GFF_FIELDS = ["seqID", "source", "compartment", "start", "end", "score", "strand", "phase", "attributes"]
COVERAGE_FIELDS = ["seqID", "start", "end", "compartment", "counts", "totalHits", "bpCovering", "compartmentLength", "coverage"]

rule extractCoverage:
    input:
        '%s/schedule_%s.json' % (out, total_buckets)
    output:
        '%s/%s/coverage/coverage_bucket_{bucket}.%s.coverage' % (out, mode, mode)
    params:
        out=Path(config['out']).resolve(),
        gff_parent=Path(config['gff_parent']).resolve(),
        extraction_parent=Path(config['extraction_parent']).resolve(),
        tempdir=config['tempdir'],
        bedtools_path=config['bedtools_path'],
        split_category=config['split_category'],
        split_collection=config['split_collection'],
    run:
        extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])
        all_compartments = [
                            ("region", None),
                            ("gene", None),
                            ("gene", "protein_coding"),
                            ("gene", "non_coding"),
                            ("exon", None),
                            ("CDS", None),
            ]
        extractions = {extract_id(file): file for file in params.extraction_parent.glob("*.csv")}
        accessions = load_bucket(wildcards.bucket)
        util_cols = ["seqID", "start", "end"]
        split_category_collection = list(map(str, params.category_collection)) + ["all"]
        coverage_table = []

        gff_cleaner = GFFCleaner(
                            tempdir=params.tempdir,
                             bedtools_path=params.bedtools_path
                            )

        for accession in accessions:
            accession = Path(accession)
            accesion_id = extract_id(accession)

            gff_file = params.gff_parent.joinpath(accession.name.replace("fna", "gff"))
            extraction_file = extractions[accession_id]

            if not gff_file.is_file():
                raise FileNotFoundError(f'Could not locate gff file for accession {accession}')

            # extraction file
            extract_df = pd.read_table(
                                       extraction_file, 
                                       usecols=["seqID", "start", "end", split_category]
                                    )
            
            accession_name = extract_name(accession)
            gff_df_merged = gff_cleaner.read(gff_file, add_exons=False)

            # Extractions = Overlapping <>
            # bedtoolsInteresect -> GroupBy -> Spacer Lengths + 'All'
            # bedtoolsCoverage -> 9 X times

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

        coverage_table = pd.concat(coverage_table, axis=0)
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

