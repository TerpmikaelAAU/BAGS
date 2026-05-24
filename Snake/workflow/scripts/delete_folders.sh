#!/usr/bin/bash -l
#SBATCH --job-name=SNAKE_WARP_DELETE
#SBATCH --output=job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abc@bio.aau.dk

# Exit on first error and if any variables are unset
set -eu

rm -r /home/bio.aau.dk/zl01hh/Aspergillus_Thea/SNAKE_WARP/data/Antismash
rm -r /home/bio.aau.dk/zl01hh/Aspergillus_Thea/SNAKE_WARP/data/GeneML_prediction