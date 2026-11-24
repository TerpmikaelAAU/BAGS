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

        # NCBI's way for big downloads: a small "dehydrated" package (reports and the list
        # of genome files in fetch.txt), then rehydrate fetches the genomes in parallel.
        # fetch NAME DATASETS-ARGS... unpacks into data/tmp/extracted/NAME/. NCBI downloads
        # fail now and then (dropped connection, corrupt zip, a file not fetched), so each
        # download is tried 3 times here. Not a rule-level retries: inside a SLURM job
        # Snakemake applies those at two levels and the attempts multiply.
        fetch_once() {{
            name=$1
            shift
            pkg=data/tmp/extracted/$name/ncbi_dataset
            rm -rf data/tmp/$name.zip data/tmp/extracted/$name
            datasets download genome "$@" --dehydrated --no-progressbar --filename data/tmp/$name.zip \\
                && unzip -q data/tmp/$name.zip -d data/tmp/extracted/$name/ \\
                && datasets rehydrate --directory data/tmp/extracted/$name/ --max-workers 10 --no-progressbar \\
                && [ -f $pkg/fetch.txt ] || return 1
            # rehydrate's handling of failed downloads is undocumented: every file in
            # fetch.txt (URL, size, path; CRLF line ends) must be there and not empty
            for f in $(cut -f3 $pkg/fetch.txt | tr -d '\\r'); do
                [ -s $pkg/$f ] || {{ echo "Not downloaded: $f" >&2; return 1; }}
            done
        }}
        fetch() {{
            for attempt in 1 2 3; do
                fetch_once "$@" && return 0
                echo "NCBI download of $1 failed (attempt $attempt of 3)" >&2
                [ $attempt = 3 ] || sleep 60
            done
            return 1
        }}

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
            fetch taxon accession --inputfile data/tmp/reference.txt
        elif [ -n "{params.taxon}" ]; then
            # GenBank only: RefSeq (GCF_) copies of GenBank assemblies would be analysed twice
            fetch taxon taxon "{params.taxon}" --assembly-source GenBank
        fi

        if [ -n "{params.accessions}" ]; then
            # Listed accessions are downloaded as given, GCF_ included; a versionless one gets its latest version
            printf '%s\\n' {params.accessions} > data/tmp/accessions.txt
            fetch accessions accession --inputfile data/tmp/accessions.txt --assembly-version latest
        fi

        # Flatten the directory structure. mv, not cp: same disk, so no second copy of
        # every genome (tens of GB for a big genus)
        mkdir -p {output.genome_dir}
        find data/tmp/extracted/ -name '*.fna' -exec mv -f {{}} {output.genome_dir}/ \\;

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
