#!/usr/bin/env python3
"""Whole-genome synteny circos plot between two assemblies (MUMmer).

The reference contigs fill the left half of the circle (clockwise, largest
first) and the query contigs the right half (counter-clockwise, largest
first), so both genomes read top to bottom. Contigs are numbered by size.

The alignment uses pyGenomeViz's MUMmer wrapper, which runs

    nucmer --mum ref.fna query.fna
    delta-filter -1            (1-to-1 alignments only)
    show-coords -H -T -r -k

and keeps the MUMmer output files in --workdir. Grey links = same
orientation, red links = inverted.
"""

import argparse
import os

from pycirclize import Circos
from pygenomeviz.align import MUMmer
from pygenomeviz.parser import Fasta


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ref", required=True, help="Reference genome FASTA")
    ap.add_argument("--query", required=True, help="Query genome FASTA")
    ap.add_argument("--ref-label", required=True, help="Reference name shown on the plot")
    ap.add_argument("--query-label", required=True, help="Query name shown on the plot")
    ap.add_argument("--out", required=True, help="Output image")
    ap.add_argument("--workdir", required=True, help="Directory for the MUMmer result files")
    ap.add_argument("--ticks-interval", type=int, default=1_000_000)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    ref = Fasta(args.ref)
    query = Fasta(args.query)

    ref_sizes = dict(sorted(ref.get_seqid2size().items(), key=lambda x: x[1], reverse=True))
    query_sizes = dict(sorted(query.get_seqid2size().items(), key=lambda x: x[1], reverse=True))

    circos = Circos(
        sectors={**ref_sizes, **dict(reversed(list(query_sizes.items())))},
        start=-358,
        end=2,
        space=4,
        sector2clockwise={seqid: False for seqid in query_sizes},
    )
    circos.text(
        f"{args.ref_label}\n({ref.full_genome_length / 1e6:.1f} Mb)", r=150, deg=35, size=10
    )
    circos.text(
        f"{args.query_label}\n({query.full_genome_length / 1e6:.1f} Mb)", r=150, deg=-35, size=10
    )

    # Number contigs 1..n by size in each genome.
    number = {name: i for i, name in enumerate(ref_sizes, start=1)}
    number.update({name: i for i, name in enumerate(query_sizes, start=1)})

    for sector in circos.sectors:
        track = sector.add_track((99.8, 100))
        track.axis(fc="black")
        if sector.size >= args.ticks_interval:
            track.xticks_by_interval(
                args.ticks_interval,
                label_formatter=lambda v: f"{v / 1e6:.1f} Mb",
                label_orientation="vertical",
            )
        sector.text(str(number[sector.name]), size=8, r=120, orientation="horizontal")

    os.makedirs(args.workdir, exist_ok=True)
    for ac in MUMmer([query, ref], outdir=args.workdir).run():
        color = "red" if ac.is_inverted else "grey"
        circos.link(
            (ac.query_name, ac.query_start, ac.query_end),
            (ac.ref_name, ac.ref_start, ac.ref_end),
            color=color,
            r1=98,
            r2=98,
        )

    fig = circos.plotfig()
    fig.savefig(args.out, dpi=args.dpi)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
