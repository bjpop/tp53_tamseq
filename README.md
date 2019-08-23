# Set up

```
conda env create --name tp53_tamseq --file environment.yaml
conda activate tp53_tamseq 
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
