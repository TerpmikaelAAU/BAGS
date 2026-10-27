rule Antismash_database:
    output:
        dir = directory("data/databases/antismashdatabase"),
    threads:
        12
    resources:
        mem_mb=resources["antismash_database"]["mem_mb"],
        runtime=resources["antismash_database"]["runtime"]
    container:
        ANTISMASH_CONTAINER
    shell:
        """
        download-antismash-databases --database-dir data/databases/antismashdatabase
        """