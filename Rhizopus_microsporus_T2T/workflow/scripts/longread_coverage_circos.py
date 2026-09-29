#!/usr/bin/env python3
"""Circos plot of long-read support for a group of chromosomes.

Tracks, outside in:
  outer ring  chromosome axis (50 kb minor ticks), telomeres outlined in orange
  93-97       CDS (teal)
  89-92       tRNA (black)
  45-88       individual ultra-long read alignments (sky blue), packed like a
              genome browser; stretches with no ultra-long alignment are shaded red
  30-44       read depth of the full read set in --window-size bp windows (green)

Inputs are minimap2 (-x map-ont) PAF files; every alignment line counts,
as in the published figures. An ultra-long alignment is drawn when it spans
at least --min-read-length bp of the chromosome. Depth is the number of
alignments touching each window, scaled to the highest window in the group.
"""

import argparse

import numpy as np
from circos_common import TELOMERE_COLOR, TRNA_COLOR, draw_telomeres, read_telomere_windows
from matplotlib.patches import Patch
from pycirclize import Circos
from pygenomeviz.parser import Genbank

CDS_TRACK = (93, 97)
TRNA_TRACK = (89, 92)
READ_TRACK = (45, 88)
COVERAGE_TRACK = (30, 44)
CDS_COLOR = "teal"
READ_COLOR = "skyblue"
DEPTH_COLOR = "green"
UNCOVERED_COLOR = "#F7AAA3"
READ_PACKING_GAP = 1500  # bp kept free between reads placed on the same row


def iter_paf_targets(paf, targets):
    """Yield (target, start, end) for PAF lines hitting one of `targets`."""
    with open(paf) as fh:
        for line in fh:
            parts = line.split("\t", 9)
            if len(parts) < 9 or parts[5] not in targets:
                continue
            yield parts[5], int(parts[7]), int(parts[8])


def ultralong_alignments(paf, chroms, min_length):
    reads = {c: [] for c in chroms}
    for chrom, start, end in iter_paf_targets(paf, reads):
        if end - start >= min_length:
            reads[chrom].append((start, end))
    return reads


def binned_depth(paf, sizes, window):
    depth = {c: np.zeros(size // window + 1) for c, size in sizes.items()}
    for chrom, start, end in iter_paf_targets(paf, depth):
        depth[chrom][start // window : end // window + 1] += 1
    return depth


def pack_reads(reads, gap=READ_PACKING_GAP):
    """Greedy row assignment: each read goes on the first row it fits on. Returns [(start, end, row)], n_rows."""
    row_ends, packed = [], []
    for start, end in sorted(reads):
        for row, row_end in enumerate(row_ends):
            if start > row_end + gap:
                row_ends[row] = end
                packed.append((start, end, row))
                break
        else:
            row_ends.append(end)
            packed.append((start, end, len(row_ends) - 1))
    return packed, len(row_ends)


def uncovered_intervals(reads, size):
    gaps, covered_to = [], 0
    for start, end in sorted(reads):
        if start > covered_to:
            gaps.append((covered_to, start))
        covered_to = max(covered_to, end)
    if covered_to < size:
        gaps.append((covered_to, size))
    return gaps


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--gbk", required=True, help="Annotated GenBank of the assembly the reads were mapped to"
    )
    ap.add_argument("--telomere", required=True, help="tidk search *_telomeric_repeat_windows.tsv")
    ap.add_argument("--ultralong-paf", required=True, help="minimap2 PAF of the ultra-long reads")
    ap.add_argument("--depth-paf", required=True, help="minimap2 PAF of the full read set")
    ap.add_argument("--chromosomes", nargs="+", required=True, help="Chromosomes to draw, in order")
    ap.add_argument("--title", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-read-length", type=int, default=50_000)
    ap.add_argument("--window-size", type=int, default=5_000)
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args()

    gbk = Genbank(args.gbk)
    all_sizes = gbk.get_seqid2size()
    sizes = {c: all_sizes[c] for c in args.chromosomes}
    features = gbk.get_seqid2features(feature_type=None)

    reads = ultralong_alignments(args.ultralong_paf, sizes, args.min_read_length)
    depth = binned_depth(args.depth_paf, sizes, args.window_size)
    max_depth = max(float(d.max()) for d in depth.values()) or 1

    circos = Circos(sectors=sizes, space=6)
    circos.text(args.title, size=14, r=115, weight="bold")

    for sector in circos.sectors:
        outer = sector.add_track((98, 100), name="outer")
        outer.axis(fc="lightgrey", ec="black")
        minor = np.arange(0, sector.size, 50_000)
        outer.xticks(minor, labels=[""] * len(minor), tick_length=1.0)
        outer.xticks(
            [0, sector.size],
            ["0 Mb", f"{sector.size / 10**6:.1f} Mb"],
            label_size=6,
            label_orientation="vertical",
            tick_length=2.5,
        )
        sector.text(sector.name, r=102, size=8, orientation="vertical")

        cds_track = sector.add_track(CDS_TRACK, r_pad_ratio=0.1)
        trna_track = sector.add_track(TRNA_TRACK, r_pad_ratio=0.1)
        read_track = sector.add_track(READ_TRACK, r_pad_ratio=0.0)
        read_track.axis(fc="#FAFAFA", ec="lightgrey", lw=0.5)
        depth_track = sector.add_track(COVERAGE_TRACK, r_pad_ratio=0.0)
        depth_track.axis(fc="#F9F9F9", ec="lightgrey", lw=0.5)

        for feat in features.get(sector.name, []):
            if feat.type == "CDS":
                cds_track.genomic_features(feat, fc=CDS_COLOR, lw=0)
            elif feat.type == "tRNA":
                trna_track.genomic_features(feat, fc=TRNA_COLOR, lw=0.1)

        chrom_reads = [(max(0, s), min(e, sector.size)) for s, e in reads[sector.name]]
        for start, end in uncovered_intervals(chrom_reads, sector.size):
            if end - start > 10:
                read_track.rect(
                    start, end, r_lim=READ_TRACK, fc=UNCOVERED_COLOR, ec="none", alpha=0.6
                )

        packed, n_rows = pack_reads(chrom_reads)
        if n_rows:
            row_height = (READ_TRACK[1] - READ_TRACK[0]) / n_rows
            for start, end, row in packed:
                r_top = READ_TRACK[1] - row * row_height
                read_track.rect(
                    start,
                    end,
                    r_lim=(r_top - 0.82 * row_height, r_top),
                    fc=READ_COLOR,
                    ec="none",
                    alpha=0.9,
                )

        bins = depth[sector.name]
        x = np.minimum(np.arange(len(bins)) * args.window_size, sector.size)
        depth_track.fill_between(
            x, bins, vmin=0, vmax=max_depth, fc=DEPTH_COLOR, ec=DEPTH_COLOR, lw=0.4, alpha=0.85
        )
        if max_depth > 110:
            ticks, labels = [0, 100, max_depth], ["0x", "100x", f"{int(max_depth)}x"]
        else:
            ticks, labels = [0, max_depth], ["0x", f"{int(max_depth)}x"]
        depth_track.yticks(ticks, labels=labels, vmax=max_depth, side="left", label_size=2)

    draw_telomeres(circos, read_telomere_windows(args.telomere))

    fig = circos.plotfig()
    handles = [
        Patch(color=CDS_COLOR, label="CDS"),
        Patch(color=TRNA_COLOR, label="tRNA"),
        Patch(color=TELOMERE_COLOR, label="Telomere"),
        Patch(color=UNCOVERED_COLOR, label="No ultra-long read coverage"),
        Patch(color=READ_COLOR, label="Ultra-long Reads"),
        Patch(color=DEPTH_COLOR, label="Full Dataset Read Depth"),
    ]
    circos.ax.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.5, 0.5),
        frameon=False,
        fontsize=5.5,
        labelspacing=0.25,
        handlelength=1.2,
        handletextpad=0.4,
    )
    fig.savefig(args.out, dpi=args.dpi)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
