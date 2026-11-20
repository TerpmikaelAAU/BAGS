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
if REFERENCE_ONLY and not requested_taxon():
    raise ValueError("config: Reference_only needs a Taxon")


checkpoint download_genomes:
    output:
        genome_dir=directory("data/genomes"),
    params:
        taxon=requested_taxon(),
        accessions=" ".join(requested_accessions()),
        reference_only=REFERENCE_ONLY,
    container:
        NCBI_DATASETS_CONTAINER
    resources:
        mem_mb=resources["NCBI_datasets"]["mem_mb"],
        runtime=resources["NCBI_datasets"]["runtime"],
    shell:
        """
        rm -rf data/tmp
        mkdir -p data/tmp/extracted

        if [ -n "{params.taxon}" ] && [ "{params.reference_only}" = True ]; then
            # NCBI's reference genomes (at most one per species). Some are flagged only on
            # the RefSeq (GCF_) copy, so take each one's GenBank (GCA_) twin. A twin can be an
            # older version of a GCA_ that is itself listed, so keep the newest per assembly.
            datasets summary genome taxon "{params.taxon}" --reference --assembly-source all --as-json-lines \\
                | awk '{{ if (match($0, /"accession": *"GCA_[0-9]+\\.[0-9]+"/) || match($0, /"pairedAccession": *"GCA_[0-9]+\\.[0-9]+"/)) {{ s = substr($0, RSTART, RLENGTH); sub(/.*"GCA_/, "GCA_", s); sub(/"$/, "", s); print s }} }}' \\
                | sort -t. -k1,1 -k2,2nr | awk -F. '!seen[$1]++' > data/tmp/reference.txt
            if [ ! -s data/tmp/reference.txt ]; then
                echo "No NCBI reference genomes for taxon {params.taxon}" >&2
                exit 1
            fi
            datasets download genome accession --inputfile data/tmp/reference.txt --filename data/tmp/taxon.zip
            unzip -q data/tmp/taxon.zip -d data/tmp/extracted/taxon/
        elif [ -n "{params.taxon}" ]; then
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
        # The reference list is checked too, in case a GenBank twin was not returned.
        REFERENCE=""
        if [ -s data/tmp/reference.txt ]; then REFERENCE=$(cat data/tmp/reference.txt); fi
        MISSING=""
        for acc in {params.accessions} $REFERENCE; do
            grep -qF -e "\\"accession\\":\\"$acc\\"" -e "\\"accession\\":\\"$acc." {output.genome_dir}/assembly_data_report.jsonl \\
                || MISSING="$MISSING $acc"
        done
        if [ -n "$MISSING" ]; then
            echo "Accessions not found on NCBI:$MISSING" >&2
            exit 1
        fi

        rm -rf data/tmp
        """
