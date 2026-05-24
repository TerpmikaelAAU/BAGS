#!/usr/bin/bash -l
#SBATCH --job-name=SNAKE_WARP_FIND
#SBATCH --output=job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=general
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abc@bio.aau.dk

# Exit on first error and if any variables are unset
set -eu
conda activate stuff
python3 /home/bio.aau.dk/zl01hh/Aspergillus_Thea/SNAKE_WARP/workflow/scripts/Find.py