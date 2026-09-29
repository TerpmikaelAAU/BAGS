PHYLO = config["phylogeny"]
PHYLO_DIR = f"{RESULTS}/phylogeny"


rule ufcg_profile:
    input:
        expand(f"{RESULTS}/genomes/ncbi/{{genome}}.fna", genome=PHYLO["genomes"]),
    output:
        directory(f"{PHYLO_DIR}/profile"),
    log:
        f"{RESULTS}/logs/ufcg_profile.log",
    conda:
        conda_env("ufcg")
    threads: PHYLO["threads"]
    resources:
        mem_mb=30000,
        runtime=2 * 24 * 60,
    params:
        genome_dir=f"{PHYLO_DIR}/input_genomes",
        tmp_dir=f"{PHYLO_DIR}/tmp",
    shell:
        """
        rm -rf {params.genome_dir} {params.tmp_dir}
        mkdir -p {params.genome_dir} {params.tmp_dir}
        for f in {input}; do ln -s "$(realpath "$f")" {params.genome_dir}/; done
        ufcg profile -i {params.genome_dir} -o {output} -w {params.tmp_dir} \
            -t {threads} --nocolor -q >{log} 2>&1
        rm -rf {params.tmp_dir}
        """


rule ufcg_tree:
    input:
        rules.ufcg_profile.output,
    output:
        directory(f"{PHYLO_DIR}/tree"),
    log:
        f"{RESULTS}/logs/ufcg_tree.log",
    conda:
        conda_env("ufcg")
    threads: PHYLO["threads"]
    resources:
        mem_mb=30000,
        runtime=24 * 60,
    params:
        multicopy="-c" if PHYLO.get("align_multicopy") else "",
    shell:
        "ufcg tree -i {input} -o {output} -t {threads} {params.multicopy} --nocolor > {log} 2>&1"


rule plot_ufcg_tree:
    input:
        rules.ufcg_tree.output,
    output:
        f"{RESULTS}/figures/ufcg_tree.png",
    log:
        f"{RESULTS}/logs/plot_ufcg_tree.log",
    conda:
        conda_env("circos")
    resources:
        mem_mb=2000,
        runtime=15,
    params:
        outgroup=PHYLO["outgroup"],
        highlight=" ".join(PHYLO.get("highlight", [])),
    shell:
        """
        tree=$(find {input} -name 'concatenated_gsi_*.nwk' | head -n 1)
        n_genes=$(basename "$tree" .nwk | sed 's/.*_//')
        python {SCRIPTS}/plot_ufcg_tree.py --tree "$tree" --outgroup {params.outgroup} \
            --highlight {params.highlight} --n-genes "$n_genes" --out {output} >{log} 2>&1
        """
