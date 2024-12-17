#!/bin/bash

j=${1:-1}
latency=${2:-30}

if [[ ! -n "$SSH_CONNECTION" ]];
then
  echo "Local environment detected."
  echo "Initializing bioinformatics genomic compartment enrichment analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
	snakemake --snakefile bedsnake_coverage.smk \
            --configfile config_coverage/config.yaml \
            --rerun-incomplete \
            --reason \
            --use-conda \
            --scheduler greedy \
            --keep-going \
            --latency-wait ${latency} \
            --cores $j
else
  echo "SSH Connection detected."
  echo "Initializing bioinformatics genomic compartment enrichment analysis. (mode ${mode}; cores ${j}) [SERVER]. Authored by Nikol Chantzi <3."
	snakemake --snakefile bedsnake_coverage.smk \
            --configfile config_coverage/config_server.yaml \
            --rerun-incomplete \
            --reason \
            --scheduler greedy \
            --keep-going \
            --cores $j \
            --latency-wait ${latency} \
            --cluster-config config_coverage/cluster_enrichment.yaml \
            --cluster "sbatch -p {cluster.partition} \
                --account={cluster.account} \
                -t {cluster.time} \
                --mem={cluster.mem} \
                --nodes={cluster.nodes} \
                --exclusive \
                -c {cluster.ncpus} \
                -J {cluster.jobName} \
                -o MindiCoverageJobDetails/{cluster.jobName}-%x-%j.out \
                -e MindiCoverageJobDetails/{cluster.jobName}-%x-%j.err"
            # --cluster "sbatch -p {cluster.partition} --account={cluster.account} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
fi