#!/usr/bin/env python3
"""Write the sequences of a GenBank file as FASTA, optionally leaving out some contigs."""

import argparse

from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gbk", required=True, help="Input GenBank file")
    ap.add_argument("--fasta", required=True, help="Output FASTA file")
    ap.add_argument("--exclude", nargs="*", default=[], help="Contig IDs to leave out")
    args = ap.parse_args()

    exclude = set(args.exclude)
    records = (r for r in SeqIO.parse(args.gbk, "genbank") if r.id not in exclude)
    n = SeqIO.write(records, args.fasta, "fasta")
    print(f"Wrote {n} records to {args.fasta}")


if __name__ == "__main__":
    main()
