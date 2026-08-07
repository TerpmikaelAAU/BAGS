rule busco_database:
    output:
        # BUSCO specifically looks for a 'lineages' subfolder inside the download path
        db_dir=directory("data/databases/busco_downloads/lineages/fungi_odb12")
    conda:
        "../envs/BUSCO.yml"
    threads: 1
    resources:
        mem_mb=4000,
        time="01:00:00"
    shell:
        """
        # Download the database to a centralized location
        busco --download fungi_odb12 --download_path data/databases/busco_downloads
        """