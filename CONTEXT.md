# BAGS

BAGS predicts genes, biosynthetic gene clusters and BUSCO completeness for every NCBI assembly of a taxon, and consolidates the results for comparison across species.

## Language

**Assembly**:
One NCBI genome assembly, identified by its accession (e.g. `GCA_000149305.1`) and named by its full file stem (e.g. `GCA_000149305.1_RO3_genomic`). The unit of one row and one data point.
_Avoid_: Genome, sample

**Strain**:
The biological isolate an assembly was sequenced from. One strain can have several assemblies.
_Avoid_: Isolate, genome

**Twin**:
A RefSeq (`GCF_`) copy of a GenBank (`GCA_`) assembly with identical sequence. Excluded by default, because counting it would double one assembly.
_Avoid_: Duplicate (ambiguous with BUSCO duplicated)

**NCBI metadata snapshot**:
The NCBI assembly report saved at download time, which fixes species, strain and quality flags for a run.

**Atypical assembly**:
An assembly NCBI flags as atypical (e.g. genome length too large, contaminated).

### BUSCO

**BUSCO category**:
One of Complete (C), Single-copy (S), Duplicated (D), Fragmented (F) or Missing (M), as a percentage of the lineage's BUSCOs. C = S + D.

### Biosynthetic gene clusters

**BGC region**:
One region reported by antiSMASH; the unit counted per assembly.
_Avoid_: Cluster, BGC (when a count is meant)

**Product class**:
The antiSMASH product label of a region (e.g. `NRPS`, `T1PKS`, `terpene`).

**Hybrid region**:
A BGC region with more than one product class. It counts once in each of its classes.
