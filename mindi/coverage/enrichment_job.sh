#!/bin/bash

#SBATCH --mem=6GB
#SBATCH --time=72:00:00
#SBATCH --partition=sla-prio
#SBATCH --account=izg5139_sc
#SBATCH --mail-type=END
#SBATCH --mail-user=nmc6088@psu.edu
#SBATCH --job-name=Mindi-Enrichment
#SBATCH --comment="Mindi Enrichment Pipeline"
#SBATCH --output=MindiEnrichment/mindi-%x-%j.out
#SBATCH --error=MindiEnrichment/mindi-%x-%j.err

echo "What a lovely day!"
echo "Current Time: $(date '+%Y-%m-%d %H:%M:%S')"
hostname
lscpu
free -h
export PYTHONBUFFERED=1
bash enrichment_sub.sh $1
