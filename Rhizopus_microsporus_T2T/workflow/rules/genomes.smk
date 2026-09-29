rule download_ncbi_genome:
    output:
        f"{RESULTS}/genomes/ncbi/{{genome}}.fna",
    log:
        f"{RESULTS}/logs/download_ncbi_genome/{{genome}}.log",
    conda:
        conda_env("ncbi")
    resources:
        mem_mb=2000,
        runtime=60,
    params:
        accession=lambda w: NCBI[w.genome],
    shell:
        """
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        datasets download genome accession {params.accession} --include genome \
            --filename "$tmp/genome.zip" >{log} 2>&1
        unzip -p "$tmp/genome.zip" 'ncbi_dataset/data/*/*.fna' >{output}
        grep -c '>' {output} | sed 's/^/sequences: /' >>{log}
        """


rule annotated_genome_fasta:
    input:
        gbk=lambda w: ANNOTATED[w.genome]["gbk"],
    output:
        f"{RESULTS}/genomes/annotated/{{genome}}.fasta",
    log:
        f"{RESULTS}/logs/annotated_genome_fasta/{{genome}}.log",
    conda:
        conda_env("circos")
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "python {SCRIPTS}/gbk_to_fasta.py --gbk {input.gbk} --fasta {output} > {log} 2>&1"
