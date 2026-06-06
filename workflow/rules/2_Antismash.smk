configfile: "config/config.yaml"

rule Antismash:
    input:
        a = "data/GeneML_prediction/{samples}.gff",
        b = "data/genomes/{samples}.fna",
    output:
        directory("data/Antismash/{samples}"),
    threads: 12,
    resources:
        mem_mb=resources["antismash"]["mem_mb"],
        runtime=resources["antismash"]["time"]
    conda:
        "../envs/Antismash.yml"
    shell:
        """
        antismash \
             -c {threads} -v \
             --databases "data/databases/antismashdatabase" \
             --genefinding-tool none --cc-mibig --cb-general \
             --output-dir "data/Antismash/{wildcards.samples}" \
             -t fungi \
             --genefinding-gff {input.a} {input.b}

        """