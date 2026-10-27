#!/usr/bin/bash -l
#SBATCH --job-name=BAGS
#SBATCH --output=SNAKE_Joblogs/job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00

# Exit on first error and if any variables are unset
set -eu
# Activate conda environment with only snakemake
conda activate BAGS
# Snakemake submits each rule as its own SLURM job via the executor plugin,
# so 1 CPU + 4G is enough here.

# Render a DAG to visualize the workflow (optional)
#snakemake --dag | dot -Tsvg > results/dag.svg

# Main workflow
snakemake --profile profile
# If the run crashed and the directory is locked, run once with --unlock
#snakemake --profile profile --unlock
