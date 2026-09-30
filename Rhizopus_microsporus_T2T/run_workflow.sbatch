#!/usr/bin/env bash
#SBATCH --job-name=rhizopus_t2t
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/snakemake_%j.log
#
# Snakemake controller job. It only schedules: every rule is submitted as a
# separate SLURM job through profiles/slurm and runs in its Apptainer
# container. Run from the Rhizopus_microsporus_T2T/ folder:
#
#   sbatch run_workflow.sbatch                 # everything in config/config.yaml
#   sbatch run_workflow.sbatch -n              # dry run
#   sbatch run_workflow.sbatch results/figures/whole_genome_synteny.png
#
# Extra arguments are passed on to snakemake.
set -euo pipefail

SNAKEMAKE_ENV=${SNAKEMAKE_ENV:-snakemake_rhizopus}   # conda env from workflow/envs/snakemake.yaml
eval "$(conda shell.bash hook)"
conda activate "$SNAKEMAKE_ENV"

command -v apptainer >/dev/null || { echo "apptainer not found on $(hostname)" >&2; exit 1; }

snakemake --profile profiles/slurm "$@"
