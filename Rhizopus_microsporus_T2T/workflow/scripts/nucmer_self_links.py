#!/usr/bin/env python3
"""Self-align a genome with MUMmer and write the duplicated blocks as a link table.

    nucmer --maxmatch -l <minmatch> genome.fasta genome.fasta
    show-coords -H -T -r genome.delta

--maxmatch keeps every exact match, including repeats, so each copy of a
duplicated region is reported rather than only the best, unique hit. The
trivial self-hit (a contig aligned to itself at the same start) is removed,
and only blocks of at least --min-length bp (reference side) are kept.

Output columns: R S1 E1 Q S2 E2 L1 (tab-separated, 1-based coordinates as
reported by show-coords; S2 > E2 means the second copy is inverted).
"""

import argparse
import os
import subprocess
import sys

from circos_common import LINK_COLUMNS, log


def run_nucmer(fasta, prefix, minmatch):
    cmd = ["nucmer", "--maxmatch", f"--minmatch={minmatch}", f"--prefix={prefix}", fasta, fasta]
    log("+ " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return f"{prefix}.delta"


def iter_links(delta):
    cmd = ["show-coords", "-H", "-T", "-r", delta]
    log("+ " + " ".join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        # S1 E1 S2 E2 LEN1 LEN2 %IDY REF QUERY
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        s1, e1, s2, e2, l1 = (int(x) for x in parts[:5])
        r_name, q_name = parts[7], parts[8]
        if r_name == q_name and s1 == s2:
            continue
        yield r_name, s1, e1, q_name, s2, e2, l1
    if proc.wait() != 0:
        sys.exit(f"show-coords failed on {delta}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True, help="Output link table (TSV)")
    ap.add_argument(
        "--minmatch", type=int, default=20, help="nucmer -l (default: 20, nucmer's default)"
    )
    ap.add_argument(
        "--min-length",
        type=int,
        default=15_000,
        help="Keep blocks whose reference-side length is at least this (default: 15000)",
    )
    ap.add_argument(
        "--keep-delta", action="store_true", help="Keep the nucmer .delta file next to --out"
    )
    args = ap.parse_args()

    prefix = os.path.splitext(args.out)[0]
    delta = run_nucmer(args.fasta, prefix, args.minmatch)

    n_total = n_kept = 0
    with open(args.out, "w") as out:
        out.write("\t".join(LINK_COLUMNS) + "\n")
        for link in iter_links(delta):
            n_total += 1
            if link[-1] >= args.min_length:
                out.write("\t".join(map(str, link)) + "\n")
                n_kept += 1
    log(f"{n_kept} of {n_total} self-alignment blocks are >= {args.min_length} bp -> {args.out}")

    if not args.keep_delta:
        os.remove(delta)


if __name__ == "__main__":
    main()
