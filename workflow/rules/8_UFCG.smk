# UFCG (universal fungal core genes): profile-pro finds the 61 core marker
# genes in each assembly's geneML proteins (MMseqs2, no AUGUSTUS), tree builds
# a maximum-likelihood phylogeny from them (MAFFT + IQ-TREE). Runs only with
# UFCG: true in config/config.yaml.
rule ufcg_profile:
    input:
        proteins="data/GeneML_prediction/{genome}.faa",
    output:
        # UFCG names the profile after the input file minus its extension
        ucg="data/UFCG/profiles/{genome}/{genome}.ucg",
    params:
        out_dir="data/UFCG/profiles/{genome}",
        tmp_dir="data/UFCG/tmp/{genome}",
        faa="data/UFCG/tmp/{genome}/{genome}.faa",
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
        # -w drops UFCG's random temp-file prefix, so each assembly needs its own dir
        rm -rf {params.tmp_dir}
        mkdir -p {params.out_dir} {params.tmp_dir}/work
        # profile-pro (1.0.6) never reads the last record of its input, so end
        # the copy with an empty one. geneML runs with --max-transcripts 1, so
        # isoforms can't make a marker look duplicated (ufcg tree drops those).
        cat {input.proteins} > {params.faa}
        echo ">end" >> {params.faa}
        ufcg profile-pro -i {params.faa} -o {params.out_dir} -w {params.tmp_dir}/work \
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
        for f in {input}; do cp "$f" data/UFCG/tree_input/; done

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
