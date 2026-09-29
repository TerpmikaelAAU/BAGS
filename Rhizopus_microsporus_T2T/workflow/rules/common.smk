from pathlib import Path

RESULTS = config["results"]
SCRIPTS = Path(workflow.basedir) / "scripts"
ANNOTATED = config["annotated_genomes"]
NCBI = config["ncbi_genomes"]
LONG_READS = config["long_reads"]
READ_GROUPS = {"_".join(group): group for group in LONG_READS["chromosome_groups"]}


wildcard_constraints:
    genome="[^/]+",
    motif="[ACGTacgt]+",


def conda_env(name):
    """A .yaml path (relative to workflow/rules/) or the name of an existing conda env."""
    return config["conda_envs"][name]


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
