# Set up

```
conda env create --name tp53_tamseq --file environment.yaml
conda activate tp53_tamseq 
```

May want to start tmux or join existing tmux session

```
# start
tmux new-session -s example_name 
# re-attach
tmux attach -t example_name
```

```
conda deactivate 
```

# Clean up and remove every computed output (WARNING dangerous)!

```
./clean.sh
```

# Dry run, find out what the pipeline will do

```
snakemake -np
```

# Run locally

```
snakemake
```

Force re-run from the start:

```
snakemake --forceall
```

# Run on the cluster

```
./run_cluster.sh
```

# Update the conda environment from environment.yaml

```
conda env update --file environment.yaml
```

# If you add new files and want to run just on those:

```
snakemake -R `snakemake --list-input-changes`
```
