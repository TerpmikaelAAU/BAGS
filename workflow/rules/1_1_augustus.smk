configfile: "config/config.yaml"

rule Augustus_prediction:
    input:
        a = "data/genomes/{samples}.fna"
    output:
        a = "data/Augustus_prediction/{samples}.gff"
    threads:
        1
    resources:
        mem_mb=resources["Augustus"]["mem_mb"],
        runtime=resources["Augustus"]["time"]
    conda:
        "../envs/augustus.yml"
    shell:
        """
        augustus \
            --strand=both \
            --gff3=on \
            --stopCodonExcludedFromCDS=false \
            --species=aspergillus_oryzae \
            {input} > {output}

        """

