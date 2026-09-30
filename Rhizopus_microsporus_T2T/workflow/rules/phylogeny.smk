PHYLO = config["phylogeny"]
PHYLO_DIR = f"{RESULTS}/phylogeny"
# UFCG looks for its database in a `config/` folder next to ufcg.jar. The jar in
# the image is read-only, so a copy runs from here with the database beside it.
UFCG_HOME = f"{PHYLO_DIR}/ufcg"
UFCG_PAYLOAD = "https://ufcg.steineggerlab.workers.dev/payload"


rule ufcg_resources:
    output:
        directory(f"{UFCG_HOME}/config"),
    log:
        f"{RESULTS}/logs/ufcg_resources.log",
    container:
        container_image("reads")
    resources:
        mem_mb=2000,
        runtime=60,
    params:
        local=PHYLO.get("ufcg_config") or "",
    shell:
        """
        rm -rf {output} && mkdir -p {output}
        if [ -n "{params.local}" ]; then
            cp -rL "{params.local}"/. {output}/ 2>{log}
        else
            for part in config core; do
                curl -sSfL --retry 5 --retry-delay 60 {UFCG_PAYLOAD}/$part.tar.gz 2>>{log} \
                    | tar -xz -C {UFCG_HOME} 2>>{log}   # unpacks into config/
            done
        fi
        test -f {output}/ppx.cfg -a -d {output}/seq -a -d {output}/model
        """


rule ufcg_profile:
    input:
        genomes=expand(f"{RESULTS}/genomes/ncbi/{{genome}}.fna", genome=PHYLO["genomes"]),
        ufcg=rules.ufcg_resources.output,
    output:
        directory(f"{PHYLO_DIR}/profile"),
    log:
        f"{RESULTS}/logs/ufcg_profile.log",
    container:
        container_image("ufcg")
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
        cp /opt/conda/share/ufcg-*/ufcg.jar {UFCG_HOME}/ufcg.jar
        for f in {input.genomes}; do ln -s "$(realpath "$f")" {params.genome_dir}/; done
        java -jar {UFCG_HOME}/ufcg.jar profile -i {params.genome_dir} -o {output} -w {params.tmp_dir} \
            -t {threads} --nocolor -q >{log} 2>&1
        rm -rf {params.tmp_dir}
        """


rule ufcg_tree:
    input:
        profile=rules.ufcg_profile.output,
        ufcg=rules.ufcg_resources.output,
    output:
        directory(f"{PHYLO_DIR}/tree"),
    log:
        f"{RESULTS}/logs/ufcg_tree.log",
    container:
        container_image("ufcg")
    threads: PHYLO["threads"]
    resources:
        mem_mb=30000,
        runtime=24 * 60,
    params:
        multicopy="-c" if PHYLO.get("align_multicopy") else "",
    shell:
        """
        cp /opt/conda/share/ufcg-*/ufcg.jar {UFCG_HOME}/ufcg.jar
        java -jar {UFCG_HOME}/ufcg.jar tree -i {input.profile} -o {output} -t {threads} \
            {params.multicopy} --nocolor >{log} 2>&1
        """


rule plot_ufcg_tree:
    input:
        rules.ufcg_tree.output,
    output:
        f"{RESULTS}/figures/ufcg_tree.png",
    log:
        f"{RESULTS}/logs/plot_ufcg_tree.log",
    container:
        container_image("circos")
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
