WGS = config["whole_genome_synteny"]


rule synteny_query_fasta:
    input:
        gbk=ANNOTATED[WGS["query"]]["gbk"],
    output:
        f"{RESULTS}/whole_genome_synteny/{WGS['query']}.fasta",
    log:
        f"{RESULTS}/logs/synteny_query_fasta.log",
    container:
        container_image("circos")
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        exclude=" ".join(WGS.get("query_exclude", [])),
    shell:
        """
        python {SCRIPTS}/gbk_to_fasta.py --gbk {input.gbk} --fasta {output} \
            --exclude {params.exclude} >{log} 2>&1
        """


rule whole_genome_synteny:
    input:
        ref=f"{RESULTS}/genomes/ncbi/{WGS['reference']}.fna",
        query=rules.synteny_query_fasta.output,
    output:
        f"{RESULTS}/figures/whole_genome_synteny.png",
    log:
        f"{RESULTS}/logs/whole_genome_synteny.log",
    container:
        container_image("circos")
    resources:
        mem_mb=16000,
        runtime=240,
    params:
        workdir=f"{RESULTS}/whole_genome_synteny/mummer",
        ref_label=WGS["reference_label"],
        query_label=WGS["query_label"],
        dpi=WGS["dpi"],
    shell:
        """
        python {SCRIPTS}/whole_genome_synteny.py --ref {input.ref} --query {input.query} \
            --ref-label "{params.ref_label}" --query-label "{params.query_label}" \
            --workdir {params.workdir} --dpi {params.dpi} --out {output} >{log} 2>&1
        """
