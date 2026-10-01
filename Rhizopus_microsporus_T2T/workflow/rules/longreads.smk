READS_DIR = f"{RESULTS}/long_reads"


rule download_reads:
    output:
        f"{READS_DIR}/{{run}}.fastq.gz",
    log:
        f"{RESULTS}/logs/download_reads/{{run}}.log",
    container:
        container_image("reads")
    resources:
        mem_mb=2000,
        runtime=24 * 60,
    shell:
        """
        urls=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={wildcards.run}&result=read_run&fields=fastq_ftp&format=tsv" \
            | tail -n +2 | cut -f2 | tr ';' ' ')
        echo "ENA files for {wildcards.run}: $urls" >{log}
        [ -n "$urls" ] || {{ echo "No FASTQ found on ENA for {wildcards.run}" >> {log}; exit 1; }}
        for url in $urls; do
            curl -sSfL --retry 5 --retry-delay 30 "https://$url" >>{output}.part 2>>{log}
        done
        mv {output}.part {output}
        """


def raw_reads(wildcards):
    return LONG_READS.get("fastq") or f"{READS_DIR}/{LONG_READS['run']}.fastq.gz"


rule ultralong_reads:
    input:
        raw_reads,
    output:
        f"{READS_DIR}/ultralong.fastq.gz",
    log:
        f"{RESULTS}/logs/ultralong_reads.log",
    container:
        container_image("reads")
    threads: 8
    resources:
        mem_mb=8000,
        runtime=6 * 60,
    params:
        min_length=LONG_READS["ultralong_min_length"],
    shell:
        """
        chopper --minlength {params.min_length} --threads {threads} --input {input} 2>{log} \
            | pigz -p {threads} >{output}
        """


rule depth_reads:
    input:
        raw_reads,
    output:
        f"{READS_DIR}/min{LONG_READS['depth_min_length']}.fastq.gz",
    log:
        f"{RESULTS}/logs/depth_reads.log",
    container:
        container_image("reads")
    threads: 8
    resources:
        mem_mb=8000,
        runtime=6 * 60,
    params:
        min_length=LONG_READS["depth_min_length"],
    shell:
        "seqkit seq --min-len {params.min_length} --threads {threads} {input} -o {output} > {log} 2>&1"


rule map_long_reads:
    input:
        ref=annotated_fasta(LONG_READS["genome"]),
        reads=f"{READS_DIR}/{{readset}}.fastq.gz",
    output:
        f"{READS_DIR}/{{readset}}.paf",
    log:
        f"{RESULTS}/logs/map_long_reads/{{readset}}.log",
    container:
        container_image("reads")
    threads: 32
    resources:
        mem_mb=64000,
        runtime=24 * 60,
    shell:
        "minimap2 -t {threads} -x map-ont {input.ref} {input.reads} -o {output} > {log} 2>&1"


rule longread_coverage_circos:
    input:
        gbk=ANNOTATED[LONG_READS["genome"]]["gbk"],
        telomere=telomere_windows(LONG_READS["genome"]),
        ultralong=f"{READS_DIR}/ultralong.paf",
        depth=f"{READS_DIR}/min{LONG_READS['depth_min_length']}.paf",
    output:
        f"{RESULTS}/figures/longread_coverage/{{group}}.png",
    log:
        f"{RESULTS}/logs/longread_coverage_circos/{{group}}.log",
    container:
        container_image("circos")
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        chromosomes=lambda w: " ".join(READ_GROUPS[w.group]),
        title=ANNOTATED[LONG_READS["genome"]]["title"],
        min_length=LONG_READS["ultralong_min_length"],
        window=LONG_READS["depth_window"],
        dpi=config["plot_dpi"],
    shell:
        """
        python {SCRIPTS}/longread_coverage_circos.py --gbk {input.gbk} --telomere {input.telomere} \
            --ultralong-paf {input.ultralong} --depth-paf {input.depth} \
            --chromosomes {params.chromosomes} --title "{params.title}" \
            --min-read-length {params.min_length} --window-size {params.window} \
            --dpi {params.dpi} --out {output} >{log} 2>&1
        """
