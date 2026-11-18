"""Build an assembly x gene cluster family (GCF) table from BiG-SCAPE output and
a UPGMA dendrogram of the assemblies on Jaccard distance between their GCF sets.

usage: bgc_tree.py --bigscape-dir data/bigscape --cutoff 0.3 --matrix OUT.tsv --tree OUT.nwk GENOME.fna...

BiG-SCAPE input files are named <assembly>__<region file>.gbk in
<bigscape-dir>/input (rule bigscape), so the GBK column of its clustering tables
maps each record to an assembly. Assemblies with no BGC regions are kept as
all-zero rows.
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict


def read_families(bigscape_dir, cutoff, genomes):
    """Return {assembly: {family: number of regions}} from every bin's clustering table."""
    counts = {g: defaultdict(int) for g in genomes}
    input_dir = os.path.join(bigscape_dir, "input")
    if not os.listdir(input_dir):
        return counts  # no assembly has a BGC region; BiG-SCAPE was skipped

    # BiG-SCAPE names folders with the float cutoff (c0.3, c1.0)
    pattern = os.path.join(bigscape_dir, "output", "output_files", f"*_c{cutoff}", "*", f"*_clustering_c{cutoff}.tsv")
    tables = sorted(glob.glob(pattern))
    if not tables:
        sys.exit(f"No BiG-SCAPE clustering tables match {pattern}")

    # Longest match first, in case one assembly name is a prefix of another
    by_length = sorted(genomes, key=len, reverse=True)
    seen = set()
    for table in tables:
        with open(table) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                # A record appears once per bin it is in; count each record once per family
                key = (row["Record"], row["Family"])
                if key in seen:
                    continue
                seen.add(key)
                genome = next((g for g in by_length if row["GBK"].startswith(g + "__")), None)
                if genome is None:
                    sys.exit(f"{table}: GBK {row['GBK']} matches no assembly")
                counts[genome][row["Family"]] += 1
    return counts


def jaccard(a, b):
    if not a and not b:
        return 0.0
    return 1 - len(a & b) / len(a | b)


def newick_name(name):
    return name if re.fullmatch(r"[A-Za-z0-9_.\-]+", name) else "'" + name.replace("'", "''") + "'"


def upgma(names, dist):
    """UPGMA on a full distance matrix; returns a Newick string with branch lengths."""
    # cluster id -> (newick, size, height)
    clusters = {i: (newick_name(n), 1, 0.0) for i, n in enumerate(names)}
    d = {(i, j): dist[i][j] for i in clusters for j in clusters if i < j}
    next_id = len(names)
    while len(clusters) > 1:
        (i, j), dij = min(d.items(), key=lambda kv: (kv[1], kv[0]))
        (ni, si, hi), (nj, sj, hj) = clusters.pop(i), clusters.pop(j)
        h = dij / 2
        node = f"({ni}:{max(h - hi, 0):.6f},{nj}:{max(h - hj, 0):.6f})"
        new_d = {}
        for k in clusters:
            dik = d[(min(i, k), max(i, k))]
            djk = d[(min(j, k), max(j, k))]
            new_d[k] = (dik * si + djk * sj) / (si + sj)
        d = {key: v for key, v in d.items() if i not in key and j not in key}
        for k, v in new_d.items():
            d[(k, next_id)] = v
        clusters[next_id] = (node, si + sj, h)
        next_id += 1
    return next(iter(clusters.values()))[0] + ";"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bigscape-dir", required=True)
    parser.add_argument("--cutoff", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("genomes", nargs="+", help="Assembly .fna files (their stems are the assembly names)")
    args = parser.parse_args()

    genomes = sorted(os.path.basename(g).removesuffix(".fna") for g in args.genomes)
    counts = read_families(args.bigscape_dir, float(args.cutoff), genomes)
    families = sorted({fam for c in counts.values() for fam in c})

    with open(args.matrix, "w") as f:
        f.write("\t".join(["assembly", *families]) + "\n")
        for g in genomes:
            f.write("\t".join([g, *(str(counts[g].get(fam, 0)) for fam in families)]) + "\n")

    sets = [set(counts[g]) for g in genomes]
    if len(genomes) < 2:
        tree = f"({newick_name(genomes[0])}:0.0);"
    else:
        dist = [[jaccard(a, b) for b in sets] for a in sets]
        tree = upgma(genomes, dist)
    with open(args.tree, "w") as f:
        f.write(tree + "\n")


if __name__ == "__main__":
    main()
