configfile: "config/config.yaml"

rule Antismash_database:
    output:
        dir = directory("data/databases/antismashdatabase"),
    threads:
        12
    resources:
        mem_mb=resources["antismash_database"]["mem_mb"],
        runtime=resources["antismash_database"]["time"]
    conda:
        "../envs/Antismash.yml"
    shell:
        """
        download-antismash-databases --database-dir data/databases/antismashdatabase
        """