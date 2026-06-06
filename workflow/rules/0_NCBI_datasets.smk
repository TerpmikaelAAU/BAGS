configfile: "config/config.yaml"

rule Ncbi_datasets:
    output:
        
    threads:
        12
    resources:
        mem_mb=resources["Ncbi_datasets"]["mem_mb"],
        runtime=resources["Ncbi_datasets"]["time"]
    conda:
        "../envs/NCBI_datasets.yml"
    shell:
        """
         datasets download genome taxon ""
        """