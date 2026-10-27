# Busco Antismash GeneML Snakemake (BAGS)

Predicts genes (GeneML), biosynthetic gene clusters (antiSMASH) and fungal BUSCO completeness for NCBI assemblies chosen by taxon and/or accession. It also builds a mashtree tree (`data/phylogeny/taxon_tree.nwk`) and consolidates the results per assembly. Terms are defined in [CONTEXT.md](CONTEXT.md).

The *Rhizopus microsporus* CBS 337.62 genome figures are in [Rhizopus_T2T_paper](https://github.com/TerpmikaelAAU/Rhizopus_T2T_paper).

## Install
Tools are installed by Snakemake from the conda envs in `workflow/envs/`. Snakemake itself comes from `environment.yml`:

```
conda env create -n BAGS -f environment.yml
```

## Configure and run
Choose the assemblies in `config/config.yaml` (see its comments), and the SLURM partition and limits in `profile/config.yaml`. Then run `mkdir -p SNAKE_Joblogs` and submit `slurm_submit.sbatch`. Without that folder SLURM fails with no log. Or run directly:

```
snakemake --profile profile -s Snakefile --conda-frontend conda --use-conda
```

## Consolidate and plot
These can be run on a partial or copied `data/` folder. The `busco_done` and `antismash_done` columns show what has finished.

```
python workflow/scripts/consolidate.py --data-dir data --ncbi-report data/genomes/assembly_data_report.jsonl
python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv --group-by species --exclude-atypical-warning contaminated --out busco_species
```

Use `--exclude-atypical-warning contaminated` rather than dropping all atypical assemblies, so that "genome length too large" ones are kept. For other options, see `plot_busco.py --help` (conda env: `workflow/envs/plotting.yml`).

## Caveats
- **Snakemake 7 only.** `environment.yml` pins 7.18.2. Snakemake 9.27.1 was tested and rejects both the profile (`--cluster-status` and `--cluster-cancel` were removed) and the `"HH:MM:SS"` runtimes (the `resources` dict in `Snakefile` and the default in `profile/config.yaml`).
- **Fungi only.** The BUSCO lineage is hard-coded to `fungi_odb12`.
- **Only the first `Taxon` entry is used.** Any further entries are ignored without a warning.
- **Taxon downloads are GenBank (`GCA_`) only.** Otherwise each RefSeq (`GCF_`) twin would be analysed a second time. Listed accessions are downloaded as given, `GCF_` included.
- **Versionless accessions.** When `datasets` is given several accessions in a file, a versionless one also returns report lines for its previous versions, which have no genome files. `--assembly-version latest` prevents this.
- **A listed accession that NCBI does not return stops the run**, e.g. a typo or a suppressed assembly.
- **Changing `Taxon` or the accessions re-runs everything.** `data/genomes` is downloaded again, so every genome gets a new timestamp and is re-analysed.
