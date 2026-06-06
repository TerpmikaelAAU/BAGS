import os
from snakemake.utils import min_version

min_version("7.18.2")

configfile: "config/config.yaml" 

GB = 1024

resources = {
    "antismash_database": {"mem_mb": 4 * GB, "time": "04:00:00"},
    "GeneML": {"mem_mb": 96 * GB, "time": "03:00:00"},
    "antismash": {"mem_mb": 64 * GB, "time": "03:00:00"},
    "NCBI_datasets": {"mem_mb": 4 * GB, "time": "04:00:00"},
    "busco": {"mem_mb": 8 * GB, "time": "01:00:00"},
}

# -----------------------------------------------------------------------------
# DYNAMIC EVALUATION FUNCTIONS
# -----------------------------------------------------------------------------
def aggregate_results(wildcards):
    genome_dir = checkpoints.download_genomes.get(**wildcards).output.genome_dir
    genomes, = glob_wildcards(os.path.join(genome_dir, "{genome}.fna"))
    return expand("data/summary/{genome}_summary.tsv", genome=genomes)

# NEW FUNCTION: Gathers the raw genomes for the tree
def get_all_genomes(wildcards):
    genome_dir = checkpoints.download_genomes.get(**wildcards).output.genome_dir
    genomes, = glob_wildcards(os.path.join(genome_dir, "{genome}.fna"))
    return expand("data/genomes/{genome}.fna", genome=genomes)

# -----------------------------------------------------------------------------
# TARGET RULE
# -----------------------------------------------------------------------------
rule all:
    input:
        "data/databases/antismashdatabase",
        "data/final_BAGS_report.tsv",
        "data/phylogeny/taxon_tree.nwk"  # Require the tree

# -----------------------------------------------------------------------------
# INCLUDES
# -----------------------------------------------------------------------------
include: "workflow/rules/0_Antismash_database.smk"
include: "workflow/rules/0_NCBI_datasets.smk"
include: "workflow/rules/1_geneML.smk"
include: "workflow/rules/2_Antismash.smk"
include: "workflow/rules/3_BUSCO.smk"
include: "workflow/rules/4_Results_summary.smk"
include: "workflow/rules/5_Phylogeny.smk" # Include the new rule





