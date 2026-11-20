"""Rename tree leaves from assembly names (GCA_002083745.1_Rhimi_ATCC52814_1_genomic)
to 'Organism strain (accession)', e.g. 'Rhizopus microsporus ATCC 52814 (GCA_002083745.1)',
using the NCBI assembly report that download_genomes keeps.

usage: label_trees.py --ncbi-report data/genomes/assembly_data_report.jsonl IN.nwk OUT.nwk

Leaves without an accession or report entry keep their name.
"""
import argparse
import re

from consolidate import accession_of, ncbi_fields, read_ncbi_report

# An unquoted leaf: after '(' or ',', up to its branch length or the next ',' / ')'.
# Internal-node labels (support values) follow ')' and are left alone.
LEAF = re.compile(r"(?<=[(,])([^(),:;'\[\]\s]+)(?=[:,);])")


def label(fields, accession):
    name = fields["organism_name"]
    strain = fields["strain"]
    if strain and strain not in name:
        name = f"{name} {strain}"
    return f"{name} ({accession})"


def quote(text):
    return "'" + text.replace("'", "''") + "'"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ncbi-report", required=True)
    parser.add_argument("tree_in")
    parser.add_argument("tree_out")
    args = parser.parse_args()

    labels = {r["accession"]: label(ncbi_fields(r), r["accession"]) for r in read_ncbi_report(args.ncbi_report)}

    def rename(match):
        accession = accession_of(match.group(1))
        return quote(labels[accession]) if accession in labels else match.group(1)

    with open(args.tree_in) as f:
        tree = f.read()
    with open(args.tree_out, "w") as f:
        f.write(LEAF.sub(rename, tree))


if __name__ == "__main__":
    main()
