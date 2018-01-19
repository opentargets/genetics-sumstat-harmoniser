Summary statistics harmonisation pipeline
=========================================

Pipeline to harmonise summary statistics collection

Work in progress

Todo:
- Write documentation
- Add VCF location to the `configs/config.yaml` file

### Set up environment

```
# Install dependencies into isolated environment
conda env create -n sumstat_harmoniser_pipeline --file environment.yaml

# Activate environment
source activate sumstat_harmoniser_pipeline
```

### Run pipeline

```
# Set up configuration file
nano configs/config.yaml

# Run snake (locally)
snakemake

# Run snakemake (on Sanger Farm)
bsub -q basement -J snake_master -n 2 -o output-%J.txt -e error-%J.txt scripts/bsub_wrapper.sh
```

### Input files
TODO

### Output files
TODO
