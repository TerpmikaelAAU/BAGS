# Rhizopus microsporus var. oligosporus CBS 337.62: genome figures

Code for the figures and genome analyses in

> Terp M. *et al.* A near telomere-to-telomere genome of *Rhizopus microsporus*
> var. *oligosporus* CBS 337.62. *Mycological Progress* (in preparation).

Everything in this folder is one [Snakemake](https://snakemake.readthedocs.io)
workflow. Every rule runs in a pinned Apptainer container (no conda
environments), and on a SLURM cluster every rule runs as a separate job. Run
all commands below from this folder (`cd Rhizopus_microsporus_T2T`); all
relative paths are relative to it.

The BUSCO scores and the BUSCO violin plot come from the BAGS workflow in the
[repository root](../README.md).

## Figures

| Manuscript figure | Output (`results/figures/`) | Rule | Script |
|---|---|---|---|
| BUSCO completeness (violin plot) | BAGS, see the [repository root](../README.md) | `3_BUSCO.smk` | `../workflow/scripts/plot_busco.py` |
| UFCG phylogeny with GSI values | `ufcg_tree.png` | `ufcg_resources` → `ufcg_profile` → `ufcg_tree` → `plot_ufcg_tree` | `plot_ufcg_tree.py` |
| Whole-genome synteny CBS 337.62 vs RT-3 | `whole_genome_synteny.png` | `whole_genome_synteny` | `whole_genome_synteny.py` |
| Annotation circos of CBS 337.62 | `annotation_circos/R_microsporus_CBS_337.62.png` | `nucmer_self_links`, `tidk_search` → `annotation_circos` | `annotation_circos.py` |
| Long-read support, chr_1 + chr_2 (main) and chr_3–13 (supplementary) | `longread_coverage/<chromosomes>.png` | `download_reads` → `ultralong_reads` / `depth_reads` → `map_long_reads` → `longread_coverage_circos` | `longread_coverage_circos.py` |
| Supplementary: six comparison genomes | `annotation_circos/<genome>.png` | as the annotation circos | `annotation_circos.py` |
| Supplementary: Redundans assembly | `annotation_circos/R_microsporus_CBS_337.62_Redundans.png` | as the annotation circos | `annotation_circos.py` |
| Telomere motif per genome (tables) | `results/telomere_finder/<genome>/` | `telomere_finder` | `telomere_finder.py` |

### What the plots show

* **Annotation circos.** Tracks from outside in: contig axis with telomeres
  outlined in orange; forward CDS (red); reverse CDS (blue); tRNA (black); CDS
  in antiSMASH biosynthetic gene clusters (green); CDS with CAZy annotation
  (purple). The inner links are self-alignment blocks of at least 15 kb from
  `nucmer --maxmatch`. **Grey** links are in the same orientation and
  **orange** links are inverted.
* **Telomeres.** A contig end is marked when a 10 kb `tidk search` window within
  the first or last 10 kb holds more than 5 copies of the motif on either strand.
* **Whole-genome synteny.** MUMmer `nucmer --mum`, `delta-filter -1`. Grey links
  are the same orientation and red links are inverted. Contigs are numbered by
  size in each genome.
* **Long-read support.** Individual alignments of reads ≥ 50 kb (sky blue);
  stretches with no such alignment are shaded red. The green track is read
  depth, counted as alignments of reads ≥ 10 kb per 5 kb window.
* **Phylogeny.** UFCG single-copy core genes, IQ-TREE on the concatenated
  protein alignment. Node labels are GSI values (the number of gene trees that
  support the split) and the tree is rooted on *R. stolonifer*.

## Input files

The only files you supply are the funannotate GenBank files (with the antiSMASH
and CAZy notes the gene tracks are drawn from) and, for the comparison genomes,
their soft-masked assemblies. They are not public. Put them at these paths
(set in `config/config.yaml`, `annotated_genomes`):

| Path | Needed for | Content |
|---|---|---|
| `data/annotation/Rhizopus_microsporus_var._oligosporus_337.62.gbk` | all figures except the tree | funannotate annotation of CBS 337.62 |
| `data/annotation/Rhizopus_microsporus.Redundans.gbk` | annotation circos, telomere finder | same, for the Redundans duplication-reduced assembly |
| `data/annotation/Rhizopus_arrhizus_Z10C7.gbk` + `_masked.fasta` | annotation circos, telomere finder | comparison genome |
| `data/annotation/Rhizopus_delemar_SE_PDA7.gbk` + `_masked.fasta` | annotation circos, telomere finder | comparison genome |
| `data/annotation/Rhizopus_microsporus_gzRhiMicr1_genomic.gbk` + `_masked.fasta` | annotation circos, telomere finder | comparison genome |
| `data/annotation/Rhizopus_microsporus_MLY36.gbk` + `_masked.fasta` | annotation circos, telomere finder | comparison genome |
| `data/annotation/Rhizopus_microsporus_var._oligosporus_RT-3.gbk` + `_masked.fasta` | annotation circos, telomere finder | comparison genome |
| `data/annotation/Rhizopus_stolonifer_PRFJ02.gbk` + `_masked.fasta` | annotation circos, telomere finder | comparison genome |

Symlinks are fine:

```bash
mkdir -p data/annotation
ln -s /path/to/funannotate_results/<genome>/predict/annotate_results/<genome>.gbk data/annotation/
```

The link targets must be visible inside the containers. On BioCloud `/home`,
`/projects`, `/raw_data` and `/databases` are; on other systems add
`--apptainer-args "--bind /path/to/data"`.

When the workflow starts it lists any of these files that are missing. To
leave a figure out, switch it off under `targets` in the config. `data/` is not
tracked by git.

Downloaded by the workflow (nothing to do):

| Data | Source | Saved as |
|---|---|---|
| CBS 337.62 and comparison assemblies | NCBI, accessions in `config/config.yaml` (`ncbi_genomes`) | `results/genomes/ncbi/<name>.fna` |
| Nanopore reads of CBS 337.62 | ENA run `long_reads.run` in the config, BioProject PRJNA1389883 | `results/long_reads/<run>.fastq.gz` |
| UFCG core-gene database (about 135 MB) | `https://ufcg.steineggerlab.workers.dev/payload/` (`config` + `core`) | `results/phylogeny/ufcg/config/` |

To use reads you already have, set `long_reads.fastq` to a local (gzipped)
FASTQ instead.

The UFCG server sometimes answers HTTP 429 (too many requests); then
`ufcg_resources` fails with that code in `results/logs/ufcg_resources.log`.
Wait and rerun, or set `phylogeny.ufcg_config` to an existing UFCG `config/`
folder (for example `<conda env>/share/ufcg-1.0.6-0/config` from a UFCG
install that has run before); it is copied instead of downloaded.

## Containers

Each tool group has one container image, listed in `config/config.yaml`
(`containers.images`). The images were built with
[Seqera Containers](https://seqera.io/containers/) from the package lists in
`workflow/containers/`, which also give the exact tool versions. They are
pulled once into `containers/<name>.sif` (about 2 GB in total):

| Image | Tools | Rules |
|---|---|---|
| `circos` | pyCirclize, pyGenomeViz, Biopython, MUMmer | plots, `nucmer_self_links`, GenBank → FASTA |
| `tidk` | tidk | `tidk_search`, `telomere_finder` |
| `reads` | minimap2, chopper, SeqKit, curl | long-read download, filtering and mapping, `ufcg_resources` |
| `ncbi` | NCBI datasets | `download_ncbi_genome` |
| `ufcg` | UFCG with AUGUSTUS, MMseqs2, MAFFT, IQ-TREE | `ufcg_profile`, `ufcg_tree` |

## Running

### 1. Install Snakemake

Only Snakemake itself is installed with conda:

```bash
conda env create -n snakemake_rhizopus -f workflow/envs/snakemake.yaml
```

### 2. Add the input files

See [Input files](#input-files).

### 3. Pull the container images (login node)

BioCloud cannot build `.sif` files inside SLURM jobs, so pull them on the
login node first:

```bash
conda activate snakemake_rhizopus
bash pull_containers.sh
```

### 4. Run on SLURM (AAU BioCloud)

The controller runs as a small SLURM job, and it submits every rule as its own
job through `profiles/slurm`:

```bash
sbatch run_workflow.sbatch -n     # dry run; see logs/snakemake_<jobid>.log
sbatch run_workflow.sbatch        # full run
sq                                # follow the jobs
```

Extra arguments are passed on to `snakemake`. For example, this builds only
one figure:

```bash
sbatch run_workflow.sbatch results/figures/whole_genome_synteny.png
```

The partition (`shared`) and the default resources are set in
`profiles/slurm/config.yaml`. The per-rule memory and time limits are set in
the rules; the longest job is the self-alignment, whose limit is
`nucmer_hours` in the config.

### 5. Run locally

With Apptainer installed, after step 3:

```bash
snakemake --sdm apptainer -c 16
```

### Run times

On 16 cores:

| Step | Time |
|---|---|
| Annotation circos plot, synteny plot, telomere search | minutes |
| `nucmer --maxmatch` self-alignment | 1–3 days per genome; about a week for the repeat-rich *R. stolonifer* |
| UFCG profile and tree of 9 genomes | about 20 minutes |

## Layout

```
config/config.yaml        genomes, accessions, parameters
workflow/Snakefile        entry point
workflow/rules/*.smk      rules, one file per analysis
workflow/scripts/*.py     plotting and helper scripts (each has --help)
workflow/containers/*.yaml  package lists of the container images
workflow/envs/snakemake.yaml  conda environment for Snakemake itself
profiles/slurm/           Snakemake SLURM profile
pull_containers.sh        pulls the container images (login node)
run_workflow.sbatch       SLURM controller job
data/annotation/          your input files (not in git)
containers/               pulled .sif images (not in git)
results/                  all outputs (not in git)
```

## Software

pyCirclize 1.9.0, pyGenomeViz 1.7.0, MUMmer 3.23, tidk 0.2.65,
UFCG 1.0.6 (with AUGUSTUS 3.5.0, MMseqs2 18.8cc5c, MAFFT 7.526, IQ-TREE 3.1.4), minimap2 2.30,
chopper 0.14.1, SeqKit 2.14.0, NCBI datasets 18.38.0, Snakemake 9.27.0, Apptainer.
Exact versions are listed in `workflow/containers/`.

## License

MIT, see [LICENSE](../LICENSE).
