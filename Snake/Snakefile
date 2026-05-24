import os
from snakemake.utils import min_version

min_version("7.18.2")

configfile: "config/config.yaml"

#Antismash -> most will run with 20 gb, and some need 50 gb.
resources = {
  "antismash_database": {"mem_mb": 4000, "time": "04:00:00"},
  "GeneML": {"mem_mb": 30000, "time": "00:30:00"},
  "antismash": {"mem_mb": 50000, "time": "03:00:00"},
  "Augustus": {"mem_mb": 5000, "time": "04:00:00"},
  }



include: "workflow/rules/0_Antismash_database.smk" 
include: "workflow/rules/1_Geneprediction.smk" 
include: "workflow/rules/2_Antismash.smk" 

include: "workflow/rules/1_1_augustus.smk"
include: "workflow/rules/2_1_Antismash_augustus.smk"


rule all:
  input:
    expand("data/databases/antismashdatabase"),
    expand("data/GeneML_prediction/{samples}.gff", samples=config["samples"]),
    expand("data/Augustus_prediction/{samples}.gff",samples=config["samples"]),
    expand("data/Antismash/{samples}/", samples=config["samples"]),
    expand("data/Antismash_augustus/{samples}/", samples=config["samples"]),
    