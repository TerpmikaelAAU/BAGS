import os
import sys
from pathlib import Path

from snakemake.exceptions import WorkflowError

RESULTS = config["results"]
SCRIPTS = Path(workflow.basedir) / "scripts"
ANNOTATED = config["annotated_genomes"]
NCBI = config["ncbi_genomes"]
LONG_READS = config["long_reads"]
READ_GROUPS = {"_".join(group): group for group in LONG_READS["chromosome_groups"]}


wildcard_constraints:
    genome="[^/]+",
    motif="[ACGTacgt]+",


def container_image(name):
    """The .sif file of a tool group, pulled by pull_containers.sh."""
    return f"{config['containers']['folder']}/{name}.sif"


def annotated_fasta(genome):
    """FASTA of the annotated GenBank file (same contig IDs as the plot)."""
    return f"{RESULTS}/genomes/annotated/{genome}.fasta"


def self_alignment_fasta(wildcards):
    return ANNOTATED[wildcards.genome].get("fasta") or annotated_fasta(wildcards.genome)


def telomere_windows(genome):
    motif = ANNOTATED[genome].get("telomere_motif")
    if not motif:
        return []
    return f"{RESULTS}/telomeres/{genome}/{genome}_{motif}_telomeric_repeat_windows.tsv"


def local_inputs():
    """(config key, path) of every file the enabled targets need that is not downloaded."""
    want = config["targets"]
    genomes = set()
    if want.get("annotation_circos") or want.get("telomere_finder"):
        genomes |= set(ANNOTATED)
    if want.get("whole_genome_synteny"):
        genomes.add(config["whole_genome_synteny"]["query"])
    if want.get("longread_coverage"):
        genomes.add(LONG_READS["genome"])
    inputs = []
    for genome in sorted(genomes):
        inputs.append((f"annotated_genomes.{genome}.gbk", ANNOTATED[genome]["gbk"]))
        if want.get("annotation_circos") and ANNOTATED[genome].get("fasta"):
            inputs.append((f"annotated_genomes.{genome}.fasta", ANNOTATED[genome]["fasta"]))
    if want.get("longread_coverage") and LONG_READS.get("fastq"):
        inputs.append(("long_reads.fastq", LONG_READS["fastq"]))
    if want.get("phylogeny") and config["phylogeny"].get("ufcg_config"):
        inputs.append(("phylogeny.ufcg_config", config["phylogeny"]["ufcg_config"]))
    return inputs


def report_missing_inputs():
    """List missing local input files once, before Snakemake reports the first of them."""
    missing = [(key, path) for key, path in local_inputs() if not os.path.exists(path)]
    if not missing:
        return
    lines = [f"  {path}    ({key})" for key, path in missing]
    print(
        f"\nMissing input files (paths are relative to {os.getcwd()}):\n"
        + "\n".join(lines)
        + "\nAdd them (symlinks are fine) or switch the figure off under `targets` in "
        "config/config.yaml. See README.md, section Input files.\n",
        file=sys.stderr,
    )


def check_containers():
    missing = [container_image(name) for name in config["containers"]["images"]]
    missing = [path for path in missing if not os.path.exists(path)]
    if missing:
        raise WorkflowError(
            "Missing container images: " + ", ".join(missing) + ". Run "
            "`bash pull_containers.sh` on a login node first (see README.md)."
        )


if workflow.is_main_process:
    check_containers()
    report_missing_inputs()


def all_targets():
    targets = []
    want = config["targets"]
    if want.get("annotation_circos"):
        targets += expand(f"{RESULTS}/figures/annotation_circos/{{genome}}.png", genome=ANNOTATED)
    if want.get("telomere_finder"):
        targets += expand(
            f"{RESULTS}/telomere_finder/{{genome}}/{{genome}}_candidate_scores.tsv",
            genome=ANNOTATED,
        )
    if want.get("whole_genome_synteny"):
        targets.append(f"{RESULTS}/figures/whole_genome_synteny.png")
    if want.get("longread_coverage"):
        targets += expand(f"{RESULTS}/figures/longread_coverage/{{group}}.png", group=READ_GROUPS)
    if want.get("phylogeny"):
        targets.append(f"{RESULTS}/figures/ufcg_tree.png")
    return targets
