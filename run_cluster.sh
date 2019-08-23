#!/bin/bash
# run on the cluster

set -o errexit

MAXJOBS=512
PREFIX=test

mkdir -p logs
echo "starting live run at $(date)..."
snakemake -p -j $MAXJOBS --latency-wait 30 --cluster-config cfg/cluster.yaml --stats logs/snakemake_stats.json --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.cores} --nodes={cluster.nodes} -t {cluster.time} --mem={cluster.memory} --output=logs/slurm-%j.out --error=logs/slurm-%j.out"

echo "finished at $(date)"

