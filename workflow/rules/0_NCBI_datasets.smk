checkpoint download_genomes:
    output:
        genome_dir=directory("data/genomes"),
    params:
        # Extracting the string from the list in your config
        taxon=config["Taxon"][0],
    conda:
        "../envs/NCBI_datasets.yml"
    resources:
        mem_mb=resources["NCBI_datasets"]["mem_mb"],
        runtime=resources["NCBI_datasets"]["time"],
    shell:
        """
        mkdir -p data/tmp
        # GenBank only: RefSeq (GCF_) copies of GenBank assemblies would be analysed twice
        datasets download genome taxon "{params.taxon}" --assembly-source GenBank --filename data/tmp/dataset.zip
        unzip -q data/tmp/dataset.zip -d data/tmp/extracted/
        
        # Flatten the directory structure
        mkdir -p {output.genome_dir}
        find data/tmp/extracted/ -name '*.fna' -exec cp {{}} {output.genome_dir}/ \\;

        # Keep NCBI's metadata snapshot (species, strain, quality flags) for consolidate.py
        cp data/tmp/extracted/ncbi_dataset/data/assembly_data_report.jsonl {output.genome_dir}/
        
        rm -rf data/tmp
        """