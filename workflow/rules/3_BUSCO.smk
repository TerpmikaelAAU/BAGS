rule busco:
    input:
        proteins="data/GeneML_prediction/{genome}.faa",
        db="data/databases/busco_downloads/lineages/fungi_odb12",
    output:
        summary="data/busco/{genome}/short_summary.txt",
        out_dir=directory("data/busco/{genome}"),
    conda:
        "../envs/BUSCO.yml"
    threads: 8
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        time=resources["busco"]["time"],
    shell:
        """
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