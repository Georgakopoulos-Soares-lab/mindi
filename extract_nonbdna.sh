#!/bin/bash

j=${1:-1}
mode=${2:-"IR"}

if [[ $(grep -i Microsoft /proc/version) ]];
then
	echo "Initializing bioinformatics nonBDNA extraction analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
	if [[ $mode == "g4" ]]; then
		snakemake --snakefile nonbdna_pipe.py --configfile config/config_g4hunter.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "regex" ]]; then
		snakemake --snakefile nonbdna_pipe.py --configfile config/config_RE.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "tandem" ]]; then
		snakemake --snakefile nonbdna_pipe.py --configfile config/config_STR.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "IR" ]]; then
		snakemake --snakefile nonbdna_pipe.py --configfile config/config_IR.yaml --latency-wait 5 --cores $j
	elif [[ $mode == "MR" ]]; then
		snakemake --snakefile nonbdna_pipe.py --configfile config/config_MR.yaml --latency-wait 5 --cores $j
  fi
else
  echo "Initializing bioinformatics genomic compartment coverage analysis. (mode ${mode}; cores ${j}) [SERVER]. Authored by Nikol Chantzi <3."
  if [[ $mode == "g4" ]]; then
    snakemake --snakefile nonbdna_pipe.py --configfile config/config_g4server.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
  elif [[ $mode == "regex" ]]; then
    snakemake --snakefile nonbdna_pipe.py --configfile config/config_RE.server.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
  elif [[ $mode == "tandem" ]]; then
    snakemake --snakefile nonbdna_pipe.py --configfile config/config_STR.server.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err"
  elif [[ $mode == "IR" ]]; then
    snakemake --snakefile nonbdna_pipe.py --configfile config/config_IR.server.yaml -j $j --latency-wait 45 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t  {cluster.time} --mem={cluster.mem} --account={cluster.account} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-IR-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-IR-%j.err"
  elif [[ $mode == "MR" ]]; then
    snakemake --snakefile nonbdna_pipe.py --configfile config/config_MR.server.yaml -j $j --latency-wait 5 --cluster-config config/cluster_overlap.yaml --cluster "sbatch -p {cluster.partition} -t  {cluster.time} --mem={cluster.mem} --account={cluster.account} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-MR-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-MR-%j.err"
  else
    echo "Request could not be fullfiled. Reason: Invalid mode $mode."
    exit 1
  fi
fi
