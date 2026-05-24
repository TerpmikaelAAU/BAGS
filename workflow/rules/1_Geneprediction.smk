configfile: "config/config.yaml"

rule GeneML_prediction:
    input:
        a = "data/genomes/{samples}.fna"
    output:
        a = "data/GeneML_prediction/{samples}.gff"
    threads:
        12
    resources:
        mem_mb=resources["GeneML"]["mem_mb"],
        runtime=resources["GeneML"]["time"]
    conda:
        "../envs/geneML.yml"
    shell:
        """
        geneml {input} -c $(nproc) -v --max-transcripts 1 --output {output}  

        """
