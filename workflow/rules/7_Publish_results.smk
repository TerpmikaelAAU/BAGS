# Copies the final outputs to results/<run start>/ (RUN_ID in the Snakefile).
# Local, because only the main Snakemake process knows this run's RUN_ID.
localrules: label_tree, publish_results


# Published trees name their leaves 'Organism strain (accession)' from NCBI's
# assembly report; data/phylogeny/ keeps the assembly-name versions.
rule label_tree:
    input:
        tree="data/phylogeny/{tree}.nwk",
        genome_dir="data/genomes",
    output:
        "data/phylogeny/labelled/{tree}.nwk",
    wildcard_constraints:
        tree="taxon_tree|bgc_tree|ufcg_tree",
    container:
        PYTHON_CONTAINER
    shell:
        """
        python workflow/scripts/label_trees.py \
            --ncbi-report {input.genome_dir}/assembly_data_report.jsonl \
            {input.tree} {output}
        """


rule publish_results:
    input:
        report="data/final_BAGS_report.tsv",
        bgc_families="data/phylogeny/bgc_families.tsv",
        trees=expand(
            "data/phylogeny/labelled/{tree}.nwk",
            tree=["taxon_tree", "bgc_tree"] + (["ufcg_tree"] if RUN_UFCG else []),
        ),
    output:
        directory(f"results/{RUN_ID}"),
    shell:
        """
        mkdir -p {output}
        cp {input} {output}/
        """
