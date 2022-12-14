#!/bin/bash
#SBATCH --partition=ccr
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL

set -e
module load python/3.8
module load snakemake/7.7.0

mkdir -p Reports

if [ -f Reports/snakemake.log ]; then
    modtime1=`stat -c %y Reports/snakemake.log|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
    mv Reports/snakemake.log Reports/snakemake.log.$modtime1
fi
touch Reports/snakemake.log

snakemake -s TETranscripts.Snakefile --latency-wait 10 --printshellcmds --cluster-config cluster.json --keep-going --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname}" -j 100 --rerun-incomplete --stats test.stats -T 2 >&1|tee -a Reports/snakemake.log

mkdir -p logs/
mv slurm-*out logs/
