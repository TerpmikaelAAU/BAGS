rule busco:
    input:
        proteins="data/GeneML_prediction/{genome}.faa",
        db="data/databases/busco_downloads/lineages/fungi_odb12",
    output:
        summary="data/busco/{genome}/short_summary.txt",
        out_dir=directory("data/busco/{genome}"),
    params:
        lineage_date=BUSCO_LINEAGE_DATE,
    container:
        BUSCO_CONTAINER
    threads: 8
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
    shell:
        """
        # An older run may have left a different fungi_odb12 release behind
        grep -qx "creation_date={params.lineage_date}" {input.db}/dataset.cfg || {{
            echo "{input.db} is not the {params.lineage_date} release; delete it to re-download" >&2
            exit 1
        }}

        # 1. Run BUSCO
        busco -i {input.proteins} -o {wildcards.genome} \
            --out_path data/busco \
            -l fungi_odb12 \
            -m proteins \
            -c {threads} \
            --offline \
            --download_path data/databases/busco_downloads \
            -f
            
        ## 2. Copy the dynamically named summary file to the static Snakemake output target
        cp {output.out_dir}/short_summary.*.txt {output.summary}
        """