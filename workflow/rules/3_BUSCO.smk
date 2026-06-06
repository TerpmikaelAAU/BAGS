rule busco:
    input:
        proteins="data/GeneML_prediction/{genome}.faa",
    output:
        summary="data/busco/{genome}/short_summary.txt",
        out_dir=directory("data/busco/{genome}"),
    conda:
        "../envs/BUSCO.yml"
    threads: 8
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["time"],
    shell:
        """
        busco -i {input.proteins} -o {wildcards.genome} \
            --out_path data/busco -l fungi_odb12 -m proteins \
            -c {threads}
        """