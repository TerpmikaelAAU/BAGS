rule antismash:
    input:
        fasta="data/genomes/{genome}.fna",
        db="data/databases/antismashdatabase",
        # Explicitly demand the GFF from the GeneML rule
        gff="data/GeneML_prediction/{genome}.gff",
    output:
        out_dir=directory("data/Antismash/{genome}"),
    container:
        ANTISMASH_CONTAINER
    threads: 8
    resources:
        mem_mb=resources["antismash"]["mem_mb"],
        runtime=resources["antismash"]["runtime"],
    shell:
        """
        antismash \
             -c {threads} -v \
             --databases {input.db} \
             --genefinding-tool none \
             --cc-mibig --cb-general \
             --output-dir {output.out_dir} \
             -t fungi \
             --genefinding-gff3 {input.gff} \
             {input.fasta}
        """