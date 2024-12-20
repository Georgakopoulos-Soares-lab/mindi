#!/bin/bash

j=${1:-1}
mode=${2:-"IR"}

if [[ ! -n "$SSH_CONNECTION" ]];
then
    echo "Initializing bioinformatics genomic compartment enrichment analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
    if [[ $mode == "g4" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_g4hunter.yaml --latency-wait 5 --cores $j
    elif [[ $mode == "regex" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_regex.yaml --latency-wait 5 --cores $j
    elif [[ $mode == "tandem" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_tandem.yaml --latency-wait 5 --cores $j
    elif [[ $mode == "IR" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_coverage.IR.yaml --latency-wait 5 --cores $j
    elif [[ $mode == "MR" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_coverage.MR.yaml --latency-wait 5 --cores $j
    else
	echo "Request could not be fullfiled. Reason: Invalid mode $mode."
	exit 1
    fi
else
    echo "Initializing bioinformatics genomic compartment enrichment analysis. (mode ${mode}; cores ${j}) [SERVER]. Authored by Nikol Chantzi <3."
    if [[ $mode == "g4" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_g4server.yaml -j $j --latency-wait 45 --cluster-config config_coverage/cluster_enrichment.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
    elif [[ $mode == "regex" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_regexserver.yaml -j $j --latency-wait 45 --cluster-config config_coverage/cluster_enrichment.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
    elif [[ $mode == "tandem" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_tandemserver.yaml -j $j --latency-wait 45 --cluster-config config_coverage/cluster_enrichment.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
    elif [[ $mode == "IR" ]]; then
    	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_coverage.IR.server.yaml -j $j --latency-wait 45 --cluster-config config_coverage/cluster_enrichment.yaml --cluster "sbatch -p {cluster.partition} --account={cluster.account} -t  {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOutEnrichment${mode}/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOutEnrichment${mode}/{cluster.jobName}-%j.err"
    elif [[ $mode == "MR" ]]; then
	snakemake --snakefile bedsnake_enrichment.py --configfile config_coverage/config_coverage.MR.server.yaml -j $j --latency-wait 5 --cluster-config config_coverage/cluster_enrichment.yaml --cluster "sbatch -p {cluster.partition} --account={cluster.account} -t  {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOutEnrichment${mode}/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOutEnrichment${mode}/{cluster.jobName}-%j.err"
    else
	echo "Request could not be fullfiled. Reason: Invalid mode $mode."
	exit 1
    fi
fi
