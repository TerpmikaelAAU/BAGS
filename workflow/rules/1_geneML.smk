# geneML is only on PyPI, so get_geneml installs the pinned version
# (GENEML_VERSION in the Snakefile) once into resources/geneml-<version>/
# inside the python image, and geneml runs from there in that same image.
rule get_geneml:
    output:
        geneml=protected(f"{GENEML_ENV}/bin/geneml"),
    params:
        env=GENEML_ENV,
        version=GENEML_VERSION,
    threads: 2
    log:
        "logs/get_geneml.log",
    resources:
        mem_mb=resources["get_geneml"]["mem_mb"],
        runtime=resources["get_geneml"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        """
        rm -rf {params.env}
        python -m venv {params.env} > {log} 2>&1
        {params.env}/bin/pip install --no-cache-dir "geneml=={params.version}" >> {log} 2>&1
        {output.geneml} --version >> {log} 2>&1
        """


rule geneml:
    input:
        # ancient: reinstalling geneML does not re-run every genome
        geneml=ancient(f"{GENEML_ENV}/bin/geneml"),
        fasta="data/genomes/{genome}.fna",
    output:
        gff="data/GeneML_prediction/{genome}.gff",
        proteins="data/GeneML_prediction/{genome}.faa",
    container:
        PYTHON_CONTAINER
    threads: 30
    resources:
        mem_mb=resources["GeneML"]["mem_mb"],
        runtime=resources["GeneML"]["runtime"],
        #gpus=1, GPU server is probally down
    shell:
        """
        {input.geneml} {input.fasta} --output {output.gff} -p {output.proteins} -c 1 --cpu-only
        """
