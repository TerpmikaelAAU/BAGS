rule generate_tree:
    input:
        # This calls a dynamic function we will add to the Snakefile
        genomes=get_all_genomes,
    output:
        tree="data/phylogeny/taxon_tree.nwk",
    conda:
        "../envs/mashtree.yml"
    threads: 8
    resources:
        mem_mb=8000,
        time="01:00:00",
    shell:
        """
        mkdir -p data/phylogeny
        
        # Count how many genomes were passed to the input
        GENOME_COUNT=$(echo {input.genomes} | wc -w)
        
        if [ "$GENOME_COUNT" -lt 3 ]; then
            echo "Not enough genomes to build a tree (found $GENOME_COUNT, need at least 3)."
            echo "(Insufficient_Genomes_$GENOME_COUNT:0.0);" > {output.tree}
        else
            # Run Mashtree using the allocated threads
            mashtree --numcpus {threads} {input.genomes} > {output.tree}
        fi
        """