rule tidk_search:
    input:
        lambda w: annotated_fasta(w.genome),
    output:
        f"{RESULTS}/telomeres/{{genome}}/{{genome}}_{{motif}}_telomeric_repeat_windows.tsv",
    log:
        f"{RESULTS}/logs/tidk_search/{{genome}}_{{motif}}.log",
    conda:
        conda_env("tidk")
    resources:
        mem_mb=4000,
        runtime=60,
    params:
        outdir=lambda w, output: Path(output[0]).parent,
        window=config["telomere"]["window"],
    shell:
        """
        tidk search --string {wildcards.motif} --window {params.window} \
            --output {wildcards.genome}_{wildcards.motif} --dir {params.outdir} {input} >{log} 2>&1
        """


rule telomere_finder:
    input:
        lambda w: annotated_fasta(w.genome),
    output:
        f"{RESULTS}/telomere_finder/{{genome}}/{{genome}}_candidate_scores.tsv",
    log:
        f"{RESULTS}/logs/telomere_finder/{{genome}}.log",
    conda:
        conda_env("tidk")
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        outdir=f"{RESULTS}/telomere_finder",
        window=config["telomere"]["window"],
    shell:
        """
        python {SCRIPTS}/telomere_finder.py --fasta {input} --name {wildcards.genome} \
            --outdir {params.outdir} --window {params.window} >{log} 2>&1
        """
