# BiG-SCAPE groups every antiSMASH region of every assembly into gene cluster
# families (GCFs). bgc_tree turns the assembly x GCF table into a dendrogram of
# shared biosynthetic potential (UPGMA on Jaccard distance), to compare with
# the mashtree tree. It is a similarity tree, not a phylogeny.
rule bigscape:
    input:
        antismash=get_antismash_dirs,
        db="data/databases/antismashdatabase",
    output:
        out_dir=directory("data/bigscape"),
    params:
        pfam=ANTISMASH_PFAM,
        cutoff=GCF_CUTOFF,
    container:
        BIGSCAPE_CONTAINER
    threads: 16
    log:
        "logs/bigscape.log",
    resources:
        mem_mb=resources["bigscape"]["mem_mb"],
        runtime=resources["bigscape"]["runtime"],
    shell:
        """
        # Region files are named after their contig, so prefix each with its
        # assembly to map BiG-SCAPE records back to assemblies.
        mkdir -p {output.out_dir}/input
        for d in {input.antismash}; do
            genome=$(basename "$d")
            for gbk in "$d"/*region*.gbk; do
                [ -e "$gbk" ] || continue
                cp "$gbk" "{output.out_dir}/input/${{genome}}__$(basename "$gbk")"
            done
        done

        # BiG-SCAPE refuses an empty input folder; bgc_tree then gives every
        # assembly an empty family set.
        if [ -z "$(ls -A {output.out_dir}/input)" ]; then
            echo "No antiSMASH regions in any assembly; skipping BiG-SCAPE" > {log}
            exit 0
        fi

        # -c: the default is every core on the node, not the allocation.
        # --exclude-gbk: the default "final" would drop an assembly named *final*.
        # No -m: MIBiG reference records would show up as families.
        # --no-trees: per-family alignment trees are not used and are slow.
        bigscape cluster \
            -i {output.out_dir}/input --input-mode flat \
            --include-gbk region --exclude-gbk __bags_exclude_nothing__ \
            -o {output.out_dir}/output -l bags \
            -p {params.pfam} -c {threads} \
            --gcf-cutoffs {params.cutoff} --include-singletons --no-trees \
            > {log} 2>&1
        """


rule bgc_tree:
    input:
        bigscape="data/bigscape",
        genomes=get_all_genomes,
    output:
        matrix="data/phylogeny/bgc_families.tsv",
        tree="data/phylogeny/bgc_tree.nwk",
    params:
        cutoff=GCF_CUTOFF,
    container:
        PYTHON_CONTAINER
    threads: 1
    resources:
        mem_mb=4 * GB,
        # UPGMA in pure Python is O(n^3): minutes for hundreds of assemblies
        runtime=240,
    shell:
        """
        python workflow/scripts/bgc_tree.py \
            --bigscape-dir {input.bigscape} --cutoff {params.cutoff} \
            --matrix {output.matrix} --tree {output.tree} \
            {input.genomes}
        """
