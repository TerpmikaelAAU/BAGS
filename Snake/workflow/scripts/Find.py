import json
import glob

out_file = "BGC0000357.5_hits.tsv"

with open(out_file, "w") as out:
    out.write("genome\tcontig\tregion\tlocus_tag\tproduct\tgene_id\n")

    for json_file in glob.glob(
        "/home/bio.aau.dk/zl01hh/Aspergillus_Thea/SNAKE_WARP/data/Antismash/**/*.json",
        recursive=True
    ):
        with open(json_file) as f:
            data = json.load(f)

        for record in data.get("records", []):
            if not record or "regions" not in record:
                continue

            genome = record.get("id", "NA")

            for region in record.get("regions", []):
                region_num = region.get("region_number", "NA")

                # index genes by id
                genes = {g["id"]: g for g in region.get("genes", [])}

                # knownclusterblast
                for hit in region.get("knownclusterblast", []):
                    if "BGC0000357.5" in hit.get("cluster", ""):
                        for gene_id in hit.get("genes", []):
                            g = genes.get(gene_id, {})
                            out.write(
                                f"{json_file}\t"
                                f"{genome}\t"
                                f"{region_num}\t"
                                f"{g.get('locus_tag','NA')}\t"
                                f"{g.get('product','NA')}\t"
                                f"{gene_id}\n"
                            )

                # mibig_hits (alternative location)
                for hit in region.get("mibig_hits", []):
                    if hit.get("accession") == "BGC0000357.5":
                        for gene_id in hit.get("genes", []):
                            g = genes.get(gene_id, {})
                            out.write(
                                f"{json_file}\t"
                                f"{genome}\t"
                                f"{region_num}\t"
                                f"{g.get('locus_tag','NA')}\t"
                                f"{g.get('product','NA')}\t"
                                f"{gene_id}\n"
                            )
