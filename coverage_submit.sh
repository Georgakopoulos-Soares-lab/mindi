#!/bin/bash

j=${1:-1}
mode=${2:-"IR"}

if [[ $(grep -i Microsoft /proc/version) ]];
then
	echo "Initializing bioinformatics genomic compartment coverage analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
	if [[ $mode == "g4" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config/config_g4hunter.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "regex" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config/config_regex.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "tandem" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config/config_tandem.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "IR" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config_coverage/config_coverage.IR.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "MR" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config_coverage/config_coverage.MR.yaml --latency-wait 5 --cores $j
	else
		echo "Request could not be fullfiled. Reason: Invalid mode $mode."
		exit 1
	fi
else
	echo "Initializing bioinformatics genomic compartment coverage analysis. (mode ${mode}; cores ${j}) [SERVER]. Authored by Nikol Chantzi <3."
	if [[ $mode == "g4" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config/config_g4server.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
	elif [[ $mode == "regex" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config/config_regexserver.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
	elif [[ $mode == "tandem" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config/config_tandemserver.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
	elif [[ $mode == "IR" ]]; then
    snakemake --snakefile bedsnake_coverage.py --configfile config_coverage/config_coverage.IR.server.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t  {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
	elif [[ $mode == "MR" ]]; then
		snakemake --snakefile bedsnake_coverage.py --configfile config_coverage/config_coverage.MR.server.yaml --latency-wait 5 --cluster-config config/cluster_overlap.yaml --cores $j
	else
		echo "Request could not be fullfiled. Reason: Invalid mode $mode."
		exit 1
	fi
fi
