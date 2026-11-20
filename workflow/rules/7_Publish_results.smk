# Copies the final outputs to results/<run start>/ (RUN_ID in the Snakefile).
# Local, because only the main Snakemake process knows this run's RUN_ID.
localrules: publish_results


rule publish_results:
    input:
        report="data/final_BAGS_report.tsv",
        taxon_tree="data/phylogeny/taxon_tree.nwk",
        bgc_tree="data/phylogeny/bgc_tree.nwk",
        bgc_families="data/phylogeny/bgc_families.tsv",
        ufcg_tree=["data/phylogeny/ufcg_tree.nwk"] if RUN_UFCG else [],
    output:
        directory(f"results/{RUN_ID}"),
    shell:
        """
        mkdir -p {output}
        cp {input} {output}/
        """
