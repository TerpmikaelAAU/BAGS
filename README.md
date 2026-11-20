# Busco Antismash GeneML Snakemake (BAGS)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24309696.svg)](https://doi.org/10.5281/zenodo.24309696)

Predicts genes (GeneML), biosynthetic gene clusters (antiSMASH) and fungal BUSCO completeness for NCBI assemblies chosen by taxon and/or accession. It also builds a mashtree tree, and a dendrogram of the assemblies by shared BGC families (BiG-SCAPE gene cluster families, UPGMA on Jaccard distance). Set `UFCG: true` in `config/config.yaml` to also build a [UFCG](https://github.com/steineggerlab/ufcg) core-gene tree (61 fungal marker genes found in the geneML proteins, IQ-TREE); with fewer than 4 assemblies it writes a star tree. Each run copies the final report, the trees (leaves named like `Rhizopus microsporus var. microsporus ATCC 52814 (GCA_002083745.1)`, from NCBI's metadata) and the assembly × BGC family table to `results/<run start time>/`.

The *Rhizopus microsporus* CBS 337.62 genome figures are in [Rhizopus_T2T_paper](https://github.com/TerpmikaelAAU/Rhizopus_T2T_paper).

## Install
Each tool runs in a pinned Apptainer container (set at the top of the `Snakefile`), so Apptainer must be on the compute nodes. GeneML has no container: the pipeline pip-installs it into `resources/`. Snakemake 9 and its SLURM executor plugin come from `environment.yml`:

```
conda env create -n BAGS -f environment.yml
```

An existing `BAGS` env from before the Snakemake 9 move must be removed (`conda env remove -n BAGS`) and created again.

## Configure and run
Choose the assemblies in `config/config.yaml` (see its comments), and the SLURM limits in `profile/config.yaml`. Then run `mkdir -p SNAKE_Joblogs` and submit `slurm_submit.sbatch`. Without that folder SLURM fails with no log. For gene prediction only (no antiSMASH, BUSCO or tree), submit `sbatch slurm_submit.sbatch geneml_only`. Or run everything directly:

```
snakemake --profile profile
```

## Consolidate and plot
These can be run on a partial or copied `data/` folder. The `busco_done` and `antismash_done` columns show what has finished.

```
python workflow/scripts/consolidate.py --data-dir data --ncbi-report data/genomes/assembly_data_report.jsonl
python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv --group-by species --exclude-atypical-warning contaminated --out busco_species
```

Use `--exclude-atypical-warning contaminated` rather than dropping all atypical assemblies, so that "genome length too large" ones are kept. For other options, see `plot_busco.py --help` (conda env: `workflow/envs/plotting.yml`).

## Citation
```bibtex
@software{terp_bags_2026,
  author    = {Terp, Mikael},
  title     = {{BAGS: Busco Antismash GeneML Snakemake}},
  year      = {2026},
  version   = {v1.0.2},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24309696},
  url       = {https://doi.org/10.5281/zenodo.24309696}
}
```

BibTeX for every tool, database and method BAGS uses is in [citations.bib](citations.bib).
