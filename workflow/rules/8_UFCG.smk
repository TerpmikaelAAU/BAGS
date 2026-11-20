# UFCG (universal fungal core genes): profile extracts the 61 core marker
# genes from each assembly, tree builds a maximum-likelihood phylogeny from
# them (MAFFT + IQ-TREE). Runs only with UFCG: true in config/config.yaml.
rule ufcg_profile:
    input:
        fasta="data/genomes/{genome}.fna",
    output:
        # UFCG names the profile after the input file minus its extension
        ucg="data/UFCG/profiles/{genome}/{genome}.ucg",
    params:
        out_dir="data/UFCG/profiles/{genome}",
        tmp_dir="data/UFCG/tmp/{genome}",
    container:
        UFCG_CONTAINER
    threads: 8
    log:
        "logs/ufcg_profile/{genome}.log",
    resources:
        mem_mb=resources["ufcg_profile"]["mem_mb"],
        runtime=resources["ufcg_profile"]["runtime"],
    shell:
        """
        # Own temp dir per assembly; UFCG's default is a shared /tmp
        rm -rf {params.tmp_dir}
        mkdir -p {params.out_dir} {params.tmp_dir}
        ufcg profile -i {input.fasta} -o {params.out_dir} -w {params.tmp_dir} \
            -t {threads} --nocolor > {log} 2>&1
        rm -rf {params.tmp_dir}
        """


rule ufcg_tree:
    input:
        get_ufcg_profiles,
    output:
        out_dir=directory("data/UFCG/tree"),
        tree="data/phylogeny/ufcg_tree.nwk",
    container:
        UFCG_CONTAINER
    threads: 32
    log:
        "logs/ufcg_tree.log",
    resources:
        mem_mb=resources["ufcg_tree"]["mem_mb"],
        runtime=resources["ufcg_tree"]["runtime"],
    shell:
        """
        # Only this run's profiles, so leftovers from an earlier taxon stay out
        rm -rf data/UFCG/tree_input
        mkdir -p data/UFCG/tree_input {output.out_dir}
        cp {input} data/UFCG/tree_input/

        # IQ-TREE needs at least 4 sequences
        n=$(ls data/UFCG/tree_input | wc -l)
        if [ "$n" -lt 4 ]; then
            echo "Only $n assemblies; UFCG tree needs at least 4. Writing a star tree." > {log}
            names=$(ls data/UFCG/tree_input | sed 's/\\.ucg$//' | paste -sd, -)
            echo "($names);" > {output.tree}
            exit 0
        fi

        ufcg tree -i data/UFCG/tree_input -o {output.out_dir} -n bags \
            -t {threads} --nocolor > {log} 2>&1
        # concatenated.nwk is the concatenated-alignment tree with assembly labels
        cp {output.out_dir}/concatenated.nwk {output.tree}
        """
