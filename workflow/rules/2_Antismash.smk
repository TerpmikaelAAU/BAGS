rule antismash:
    input:
        gff="data/GeneML_prediction/{genome}.gff",
        db="data/databases/antismashdatabase",
    output:
        out_dir=directory("data/Antismash/{genome}"),
    conda:
        "../envs/Antismash.yml"
    threads: 8
    resources:
        mem_mb=resources["antismash"]["mem_mb"],
        runtime=resources["antismash"]["time"],
    shell:
        """
        antismash --databases {input.db} --output-dir {output.out_dir} {input.gff}
        """