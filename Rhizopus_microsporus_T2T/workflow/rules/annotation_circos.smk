rule nucmer_self_links:
    input:
        fasta=self_alignment_fasta,
    output:
        f"{RESULTS}/self_alignment/{{genome}}_self_links.tsv",
    log:
        f"{RESULTS}/logs/nucmer_self_links/{{genome}}.log",
    conda:
        conda_env("circos")
    resources:
        mem_mb=50000,
        runtime=lambda w: 60 * ANNOTATED[w.genome].get("nucmer_hours", 48),
    params:
        minmatch=lambda w: ANNOTATED[w.genome].get("nucmer_minmatch", 20),
        min_length=config["self_links"]["min_length"],
    shell:
        """
        python {SCRIPTS}/nucmer_self_links.py --fasta {input.fasta} --out {output} \
            --minmatch {params.minmatch} --min-length {params.min_length} >{log} 2>&1
        """


rule annotation_circos:
    input:
        gbk=lambda w: ANNOTATED[w.genome]["gbk"],
        links=rules.nucmer_self_links.output,
        telomere=lambda w: telomere_windows(w.genome),
    output:
        f"{RESULTS}/figures/annotation_circos/{{genome}}.png",
    log:
        f"{RESULTS}/logs/annotation_circos/{{genome}}.log",
    conda:
        conda_env("circos")
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        title=lambda w: ANNOTATED[w.genome]["title"],
        exclude=lambda w: " ".join(ANNOTATED[w.genome].get("exclude", [])),
        telomere=lambda w, input: f"--telomere {input.telomere}" if input.telomere else "",
        min_link_length=config["self_links"]["min_length"],
        min_repeats=config["telomere"]["min_repeats"],
        dpi=config["plot_dpi"],
    shell:
        """
        python {SCRIPTS}/annotation_circos.py --gbk {input.gbk} --links {input.links} \
            {params.telomere} --title "{params.title}" --exclude {params.exclude} \
            --min-link-length {params.min_link_length} \
            --telomere-min-repeats {params.min_repeats} --dpi {params.dpi} \
            --out {output} >{log} 2>&1
        """
