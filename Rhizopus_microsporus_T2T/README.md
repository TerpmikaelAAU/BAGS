# Rhizopus microsporus var. oligosporus CBS 337.62: genome figures

Code for the figures and genome analyses in

> Terp M. *et al.* A near telomere-to-telomere genome of *Rhizopus microsporus*
> var. *oligosporus* CBS 337.62. *Mycological Progress* (in preparation).

Everything in this folder is one [Snakemake](https://snakemake.readthedocs.io)
workflow. Each rule has its own conda environment (`workflow/envs/`), and on a
SLURM cluster every rule runs as a separate job. Run all commands below from
this folder (`cd Rhizopus_microsporus_T2T`).

The BUSCO scores and the BUSCO violin plot come from the BAGS workflow in the
[repository root](../README.md).

## Figures

| Manuscript figure | Output (`results/figures/`) | Rule | Script |
|---|---|---|---|
| BUSCO completeness (violin plot) | BAGS, see the [repository root](../README.md) | `3_BUSCO.smk` | `../workflow/scripts/plot_busco.py` |
| UFCG phylogeny with GSI values | `ufcg_tree.png` | `ufcg_profile` → `ufcg_tree` → `plot_ufcg_tree` | `plot_ufcg_tree.py` |
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

## Data

Downloaded by the workflow:

| Data | Source |
|---|---|
| CBS 337.62 assembly | NCBI GCA_060230415.1 |
| Comparison assemblies | NCBI, accessions in `config/config.yaml` (`ncbi_genomes`) |
| Nanopore reads of CBS 337.62 | ENA/SRA run in `config/config.yaml` (`long_reads.run`), BioProject PRJNA1389883 |

Supplied by you: the funannotate GenBank files, including the antiSMASH and
CAZy notes the gene tracks are drawn from. They are not public; put them in
`data/annotation/` with the file names used in `config/config.yaml`:

| File | Content |
|---|---|
| `Rhizopus_microsporus_var._oligosporus_337.62.gbk` | funannotate annotation of CBS 337.62 (with antiSMASH and CAZy notes) |
| `Rhizopus_microsporus.Redundans.gbk` | same, for the Redundans duplication-reduced assembly |
| `<genome>.gbk`, `<genome>_masked.fasta` | funannotate annotation and soft-masked assembly of each comparison genome |

Symlinks are fine:

```bash
mkdir -p data/annotation
ln -s /path/to/funannotate_results/<genome>/predict/annotate_results/<genome>.gbk data/annotation/
```

`data/` is not tracked by git.

## Running

### 1. Install Snakemake

```bash
conda env create -n snakemake_rhizopus -f workflow/envs/snakemake.yaml
```

### 2. Add the annotation files

Put the annotation files in `data/annotation/` (see [Data](#data)) and
check `config/config.yaml`. To build only some of the figures, switch parts
off under `targets`.

### 3. Run on SLURM (AAU BioCloud)

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

### 4. Run locally

```bash
snakemake --use-conda -c 16
```

### Run times

On 16 cores:

| Step | Time |
|---|---|
| Annotation circos plot, synteny plot, telomere search | minutes |
| `nucmer --maxmatch` self-alignment | 1–3 days per genome; about a week for the repeat-rich *R. stolonifer* |
| UFCG profile of 9 genomes | hours |

## Layout

```
config/config.yaml        genomes, accessions, parameters
workflow/Snakefile        entry point
workflow/rules/*.smk      rules, one file per analysis
workflow/scripts/*.py     plotting and helper scripts (each has --help)
workflow/envs/*.yaml      conda environments
profiles/slurm/           Snakemake SLURM profile
run_workflow.sbatch       SLURM controller job
```

## Software

pyCirclize 1.9.0, pyGenomeViz 1.7.0, MUMmer 3.23, tidk 0.2.65,
UFCG 1.0.6 (with AUGUSTUS, MMseqs2, MAFFT, IQ-TREE), minimap2 2.30,
chopper 0.14.1, SeqKit 2.14.0, NCBI datasets 18.38.0, Snakemake 9.27.0.
Exact versions are pinned in `workflow/envs/`.

## License

MIT, see [LICENSE](../LICENSE).
