def requested_taxon():
    taxon = config.get("Taxon") or []
    if isinstance(taxon, str):
        taxon = [taxon]
    return taxon[0] if taxon else ""


def requested_accessions():
    """Accessions from config["Accessions"] (list) and config["Accessions_file"] (one per line, # comments)."""
    accessions = list(config.get("Accessions") or [])
    accessions_file = config.get("Accessions_file")
    if accessions_file:
        with open(accessions_file) as f:
            accessions += [line.split("#")[0] for line in f]
    accessions = [str(a).strip() for a in accessions if str(a).strip()]
    return list(dict.fromkeys(accessions))


if not requested_taxon() and not requested_accessions():
    raise ValueError("config: set Taxon, Accessions and/or Accessions_file")


checkpoint download_genomes:
    output:
        genome_dir=directory("data/genomes"),
    params:
        taxon=requested_taxon(),
        accessions=" ".join(requested_accessions()),
    conda:
        "../envs/NCBI_datasets.yml"
    resources:
        mem_mb=resources["NCBI_datasets"]["mem_mb"],
        runtime=resources["NCBI_datasets"]["time"],
    shell:
        """
        rm -rf data/tmp
        mkdir -p data/tmp/extracted

        if [ -n "{params.taxon}" ]; then
            # GenBank only: RefSeq (GCF_) copies of GenBank assemblies would be analysed twice
            datasets download genome taxon "{params.taxon}" --assembly-source GenBank --filename data/tmp/taxon.zip
            unzip -q data/tmp/taxon.zip -d data/tmp/extracted/taxon/
        fi

        if [ -n "{params.accessions}" ]; then
            # Listed accessions are downloaded as given, GCF_ included; a versionless one gets its latest version
            printf '%s\\n' {params.accessions} > data/tmp/accessions.txt
            datasets download genome accession --inputfile data/tmp/accessions.txt --assembly-version latest --filename data/tmp/accessions.zip
            unzip -q data/tmp/accessions.zip -d data/tmp/extracted/accessions/
        fi

        # Flatten the directory structure
        mkdir -p {output.genome_dir}
        find data/tmp/extracted/ -name '*.fna' -exec cp {{}} {output.genome_dir}/ \\;

        # Keep NCBI's metadata snapshot (species, strain, quality flags) for consolidate.py,
        # merged across the taxon and accession downloads
        find data/tmp/extracted/ -name assembly_data_report.jsonl -exec cat {{}} + \\
            | awk '!seen[$0]++' > {output.genome_dir}/assembly_data_report.jsonl

        # Fail if a listed accession was not downloaded (unknown, suppressed or mistyped).
        # A versionless accession (GCA_000149305) matches whichever version NCBI returned.
        MISSING=""
        for acc in {params.accessions}; do
            grep -qF -e "\\"accession\\":\\"$acc\\"" -e "\\"accession\\":\\"$acc." {output.genome_dir}/assembly_data_report.jsonl \\
                || MISSING="$MISSING $acc"
        done
        if [ -n "$MISSING" ]; then
            echo "Accessions not found on NCBI:$MISSING" >&2
            exit 1
        fi

        rm -rf data/tmp
        """
