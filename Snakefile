import os
from snakemake.utils import min_version

min_version("9.0")

configfile: "config/config.yaml" 

GB = 1024

# runtime is in minutes
resources = {
    "antismash_database": {"mem_mb": 4 * GB, "runtime": 240},
    "GeneML": {"mem_mb": 30 * GB, "runtime": 180},
    "antismash": {"mem_mb": 64 * GB, "runtime": 180},
    "NCBI_datasets": {"mem_mb": 4 * GB, "runtime": 240},
    "busco": {"mem_mb": 8 * GB, "runtime": 60},
    "get_geneml": {"mem_mb": 8 * GB, "runtime": 120},
}

# -----------------------------------------------------------------------------
# CONTAINERS
# Every tool rule runs in a pinned image (profile: software-deployment-method:
# apptainer). The summary rules only use the shell, so they run on the host.
# -----------------------------------------------------------------------------
NCBI_DATASETS_CONTAINER = "docker://staphb/ncbi-datasets:18.37.0"
# antiSMASH: the official image, not the biocontainer. antiSMASH writes
# pre-built data (e.g. bgc_seeds.hmm) into its own package directory, which
# is read-only in a container; the official image ships it pre-built.
ANTISMASH_CONTAINER = "docker://antismash/standalone-lite:8.0.4"
BUSCO_CONTAINER = "docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2"
MASHTREE_CONTAINER = "docker://quay.io/biocontainers/mashtree:1.4.6--pl5321h87e0c26_4"
PYTHON_CONTAINER = "docker://python:3.12-slim"

# geneML is only on PyPI (no container), so rule get_geneml installs it once
# into resources/ inside the python image, and geneml runs from there.
GENEML_VERSION = "1.3.1"
GENEML_ENV = f"resources/geneml-{GENEML_VERSION}"

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
include: "workflow/rules/0_BUSCO_database.smk"
include: "workflow/rules/0_Antismash_database.smk"
include: "workflow/rules/0_NCBI_datasets.smk"
include: "workflow/rules/1_geneML.smk"
include: "workflow/rules/2_Antismash.smk"
include: "workflow/rules/3_BUSCO.smk"
include: "workflow/rules/4_Results_summary.smk"
include: "workflow/rules/5_Phylogeny.smk" # Include the new rule





