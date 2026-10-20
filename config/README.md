# Configuration

`config.yaml` chooses which NCBI assemblies are analysed. Set any combination of the three keys; the assemblies from all of them are combined, and an assembly listed twice is downloaded once.

| Key | Value | Downloads |
|---|---|---|
| `Taxon` | `[name or taxid]`, or `[]` | Every GenBank (`GCA_`) assembly of the taxon. RefSeq copies are left out. |
| `Accessions` | List of `GCA_`/`GCF_` accessions, or `[]` | Exactly these assemblies. A versionless accession (`GCA_000149305`) gets its latest version. |
| `Accessions_file` | Path to a text file, or `""` | The same as `Accessions`, read from a file with one accession per line (`#` starts a comment). See `accessions.txt`. |

Taxon only:

```yaml
Taxon: [Ganoderma lucidum]
Accessions: []
Accessions_file: ""
```

Specific assemblies only:

```yaml
Taxon: []
Accessions:
  - GCA_000149305.1
  - GCA_000271645.1
Accessions_file: ""
```

The run stops with `Accessions not found on NCBI: ...` if a listed accession could not be downloaded, for example because it is mistyped or suppressed.
