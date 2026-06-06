rule geneml:
    input:
        fasta="data/genomes/{genome}.fna",
    output:
        gff="data/GeneML_prediction/{genome}.gff",
        proteins="data/GeneML_prediction/{genome}.faa",
    conda:
        "../envs/geneML.yml"
    threads: 1
    resources:
        mem_mb=resources["GeneML"]["mem_mb"],
        runtime=resources["GeneML"]["time"],
        #gpus=1, GPU server is probally down
    shell:
        """
        geneml {input.fasta} --output {output.gff} -p {output.proteins} -c 1 --cpu-only
        """