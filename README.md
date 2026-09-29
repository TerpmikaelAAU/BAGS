# Busco Antismash GeneML Snakemake (BAGS) 

## Requirements
All required tools are automatically installed by Snakemake using conda environments or singularity/apptainer containers, however Snakemake itself needs to be installed first. Load a software module with Snakemake, use a native install, or use the `environment.yml` file to create a conda environment for this particular project using fx `conda env create -n <snakemake_template> -f environment.yml`.

## Usage
Adjust the `config.yaml` files under both `config/` and `profiles/` accordingly, then simply run `snakemake --profile profiles/<subfolder>` or submit a SLURM job using the `slurm_submit.sbatch` example script.
The usage of this workflow is also described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).



## Consolidating and plotting results
After (or during) a run, `workflow/scripts/consolidate.py` collects the NCBI metadata, BUSCO and antiSMASH results into one row per assembly. Every assembly in the NCBI report gets a row, with `busco_done` / `antismash_done` showing what has finished, so it can be run on a partial or copied `data/` folder. Terms are defined in [CONTEXT.md](CONTEXT.md).

```
python workflow/scripts/consolidate.py --data-dir data --ncbi-report data/genomes/assembly_data_report.jsonl
```

This writes to `data/consolidated/`:
- `BAGS_assemblies.tsv`: NCBI metadata (species, strain, assembly level, atypical flag, N50, ...), BUSCO percentages and BGC region counts per product class. RefSeq copies of GenBank assemblies are marked in `twin_of`.
- `BAGS_bgc_regions.tsv`: one row per antiSMASH region with its product classes, a `bacterial_type` flag for bacterial product classes (e.g. endosymbiont contamination) and the region's GC content.

`workflow/scripts/plot_busco.py` draws violin plots of the BUSCO categories from `BAGS_assemblies.tsv`, by default with all assemblies of the taxon in one violin per category (conda env: `workflow/envs/plotting.yml`):

```
python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv --out busco_all
python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv --group-by species --exclude-atypical-warning contaminated --out busco_species
python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv --compare-unfiltered --min-busco-c 90 --out busco_compare
```

Filters: `--exclude-atypical-warning TEXT` (repeatable; drops assemblies whose NCBI atypical warning contains TEXT, e.g. `contaminated`, so "genome length too large" assemblies can be kept), `--min-busco-c`, `--species` (repeatable). The pooled plot colours each category (`--category-colors`, default `green-blue`; BUSCO's yellow and red for Fragmented and Missing). Each plot is saved as PNG, PDF and a TSV of the plotted values.
