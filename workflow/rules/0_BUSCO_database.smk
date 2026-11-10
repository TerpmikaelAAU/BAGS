rule busco_database:
    output:
        # BUSCO specifically looks for a 'lineages' subfolder inside the download path
        db_dir=directory("data/databases/busco_downloads/lineages/fungi_odb12")
    container:
        BUSCO_CONTAINER
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    params:
        date=BUSCO_LINEAGE_DATE,
        md5=BUSCO_LINEAGE_MD5,
    shell:
        """
        python workflow/scripts/get_busco_lineage.py fungi_odb12 {params.date} {params.md5} data/databases/busco_downloads
        """
