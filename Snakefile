import os
from snakemake.utils import min_version

min_version("7.18.2")

configfile: "config/config.yaml"

#Antismash -> most will run with 20 gb, and some need 50 gb.
resources = {
  "antismash_database": {"mem_mb": 4000, "time": "04:00:00"},
  "GeneML": {"mem_mb": 30000, "time": "00:30:00"},
  "antismash": {"mem_mb": 50000, "time": "03:00:00"},
  "NCBI_datasets": {"mem_mb": 4000, "time": "04:00:00"},
  
  }



include: "workflow/rules/0_Antismash_database.smk" 
include: "workflow/rules/0_NCBI_datasets.smk"
include: "workflow/rules/1_GeneML.smk" 
include: "workflow/rules/2_Antismash.smk" 



rule all:
  input:
    expand("data/databases/antismashdatabase"),
    expand("data/GeneML_prediction/{samples}.gff", samples=config["samples"]),
    expand("data/Antismash/{samples}/", samples=config["samples"]),
    expand("data/busco/{samples}/BUSCO/", samples=config["samples"])
    