rule antismash:
    input:
        fasta="data/genomes/{genome}.fna",
        db="data/databases/antismashdatabase",
        # Explicitly demand the GFF from the GeneML rule
        gff="data/GeneML_prediction/{genome}.gff",
    output:
        out_dir=directory("data/Antismash/{genome}"),
    conda:
        "../envs/Antismash.yml" # Remember to keep the ../../ path we fixed earlier!
    threads: 8
    resources:
        mem_mb=resources["antismash"]["mem_mb"],
        runtime=resources["antismash"]["time"],
    shell:
        """
        antismash \
             -c {threads} -v \
             --databases {input.db} \
             --genefinding-tool none \
             --cc-mibig --cb-general \
             --output-dir {output.out_dir} \
             -t fungi \
             --genefinding-gff {input.gff} \
             {input.fasta}
        """