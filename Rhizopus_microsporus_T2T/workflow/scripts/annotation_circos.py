#!/usr/bin/env python3
"""Circos plot of one annotated genome: gene tracks, telomeres and duplicated blocks.

Tracks, outside in:
  outer ring  contig axis, telomeres outlined in orange (tidk search windows)
  93-97       forward-strand CDS (red)
  89-93       reverse-strand CDS (blue)
  85-89       tRNA (black)
  81-85       CDS inside an antiSMASH biosynthetic gene cluster (green)
  77-81       CDS with a CAZy annotation (purple)
  <75         self-alignment blocks >= --min-link-length bp
              (grey = same orientation, orange = inverted)

Used for the main CBS 337.62 figure, the Redundans (duplication-reduced)
supplementary figure and the six comparison genomes.
"""

import argparse

from circos_common import (
    add_annotation_tracks,
    add_outer_axis,
    annotation_legend_handles,
    draw_self_links,
    draw_telomeres,
    log,
    read_links,
    read_telomere_windows,
    sort_by_size,
)
from pycirclize import Circos
from pygenomeviz.parser import Genbank


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gbk", required=True, help="Annotated GenBank file (funannotate output)")
    ap.add_argument("--links", required=True, help="Link table from nucmer_self_links.py")
    ap.add_argument("--telomere", help="tidk search *_telomeric_repeat_windows.tsv (optional)")
    ap.add_argument("--title", required=True)
    ap.add_argument("--out", required=True, help="Output image (.png/.pdf/.svg)")
    ap.add_argument("--exclude", nargs="*", default=[], help="Contig IDs to leave out of the plot")
    ap.add_argument("--min-link-length", type=int, default=15_000)
    ap.add_argument("--telomere-min-repeats", type=int, default=5)
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args()

    gbk = Genbank(args.gbk)
    sectors = sort_by_size(gbk.get_seqid2size(), args.exclude)
    seqid2features = gbk.get_seqid2features(feature_type=None)

    circos = Circos(sectors=sectors, space=2)
    circos.text(args.title, size=14, r=115, weight="bold")

    for sector in circos.sectors:
        add_outer_axis(sector)
        add_annotation_tracks(sector, seqid2features[sector.name])

    if args.telomere:
        draw_telomeres(circos, read_telomere_windows(args.telomere, args.telomere_min_repeats))
    else:
        log("No telomere table given; telomere track left empty")

    n_links = draw_self_links(circos, read_links(args.links), min_length=args.min_link_length)
    log(f"Links plotted: {n_links}")

    fig = circos.plotfig()
    circos.ax.legend(
        handles=annotation_legend_handles(),
        bbox_to_anchor=(1.07, 0.99),
        loc="upper right",
        fontsize=9,
    )
    fig.savefig(args.out, dpi=args.dpi)
    log(f"Saved {args.out}")


if __name__ == "__main__":
    main()
