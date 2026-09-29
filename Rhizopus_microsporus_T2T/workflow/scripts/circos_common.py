"""Helpers shared by the pyCirclize scripts in this workflow.

Kept in one place so every circos figure draws annotation tracks, telomeres
and self-alignment links the same way.
"""

import sys

import pandas as pd
from matplotlib.patches import Patch

# Colours used in the manuscript figures.
FWD_CDS_COLOR = "red"
REV_CDS_COLOR = "blue"
TRNA_COLOR = "black"
BGC_COLOR = "green"
CAZY_COLOR = "purple"
TELOMERE_COLOR = "orange"
DIRECT_LINK_COLOR = (0.4, 0.4, 0.4, 0.2)  # grey, 20 % alpha
INVERTED_LINK_COLOR = (1.0, 0.65, 0.0, 0.2)  # orange, 20 % alpha

LINK_COLUMNS = ["R", "S1", "E1", "Q", "S2", "E2", "L1"]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def sort_by_size(seqid2size, exclude=()):
    """Return {seqid: size} without `exclude`, largest sequence first."""
    kept = {k: v for k, v in seqid2size.items() if k not in set(exclude)}
    return dict(sorted(kept.items(), key=lambda x: x[1], reverse=True))


def add_outer_axis(sector, label_size=6, name_r=101):
    """Grey outer ring with 0 / length ticks and the contig name. Returns the track."""
    outer = sector.add_track((98, 100), name="outer")
    outer.axis(fc="lightgrey", ec="black")
    outer.xticks(
        [0, sector.size],
        ["0 Mb", f"{sector.size / 10**6:.1f} Mb"],
        label_size=label_size,
        label_orientation="vertical",
        tick_length=2,
    )
    sector.text(sector.name, r=name_r, size=8, orientation="vertical")
    return outer


def add_annotation_tracks(sector, features):
    """Forward CDS, reverse CDS, tRNA, BGC and CAZy tracks (radius 97 -> 77).

    BGC and CAZy membership is read from the funannotate `note` qualifiers
    ("antiSMASH:Cluster_..." and "CAZy:...").
    """
    f_cds = sector.add_track((93, 97), r_pad_ratio=0.1)
    r_cds = sector.add_track((89, 93), r_pad_ratio=0.1)
    trna = sector.add_track((85, 89), r_pad_ratio=0.1)
    bgc = sector.add_track((81, 85), r_pad_ratio=0.1)
    cazy = sector.add_track((77, 81), r_pad_ratio=0.1)

    for feat in features:
        if feat.type == "CDS":
            if feat.location.strand == 1:
                f_cds.genomic_features(feat, fc=FWD_CDS_COLOR, lw=0)
            else:
                r_cds.genomic_features(feat, fc=REV_CDS_COLOR, lw=0)
            note = "".join(feat.qualifiers.get("note", ""))
            if "antiSMASH:Cluster_" in note:
                bgc.rect(int(feat.location.start), int(feat.location.end), fc=BGC_COLOR, ec="none")
            if "CAZy:" in note:
                cazy.rect(
                    int(feat.location.start), int(feat.location.end), fc=CAZY_COLOR, ec="none"
                )
        elif feat.type == "tRNA":
            trna.genomic_features(feat, fc=TRNA_COLOR, lw=0.1)


def annotation_legend_handles():
    return [
        Patch(color=FWD_CDS_COLOR, label="Forward CDS"),
        Patch(color=REV_CDS_COLOR, label="Reverse CDS"),
        Patch(color=TRNA_COLOR, label="tRNA"),
        Patch(color=BGC_COLOR, label="BGC"),
        Patch(color=CAZY_COLOR, label="CAZy"),
        Patch(color=TELOMERE_COLOR, label="Telomere"),
    ]


def read_telomere_windows(tsv, min_repeats=5):
    """tidk `search` windows TSV, keeping windows with > min_repeats repeats on either strand."""
    df = pd.read_csv(tsv, sep="\t")
    keep = (df["forward_repeat_number"] > min_repeats) | (df["reverse_repeat_number"] > min_repeats)
    return df[keep]


def draw_telomeres(circos, telomere_df, window=10_000, track_name="outer"):
    """Outline the first/last `window` bp of a contig when a telomere-repeat window falls there.

    tidk labels each window by its end coordinate (the last window by the contig
    length), so a window <= `window` is the contig start and a window within
    `window` bp of the contig end is the contig end. Hits elsewhere
    (interstitial repeats) are reported but not drawn.
    """
    for sector in circos.sectors:
        hits = telomere_df[telomere_df["id"] == sector.name]["window"].astype(int)
        track = sector.get_track(track_name)
        if (hits <= window).any():
            track.rect(0, min(window, sector.size), fc="none", ec=TELOMERE_COLOR, lw=4)
        if (hits > sector.size - window).any():
            track.rect(
                max(0, sector.size - window), sector.size, fc="none", ec=TELOMERE_COLOR, lw=4
            )
        interior = hits[(hits > window) & (hits <= sector.size - window)]
        if len(interior):
            log(
                f"  {sector.name}: interstitial telomere-repeat windows not drawn: {list(interior)}"
            )


def read_links(path):
    """Self-alignment links written by nucmer_self_links.py (also reads the old *_links_cache.csv)."""
    sep = "," if str(path).endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, usecols=LINK_COLUMNS, dtype={"R": str, "Q": str})


def draw_self_links(circos, links, min_length=15_000, r=75):
    """Draw nucmer self-alignment blocks >= min_length bp between the plotted contigs.

    Grey = both copies in the same orientation, orange = inverted copy
    (show-coords reports an inverted query copy with S2 > E2).
    """
    names = {s.name for s in circos.sectors}
    links = links[(links["L1"] >= min_length) & links["R"].isin(names) & links["Q"].isin(names)]
    for row in links.itertuples(index=False):
        color = INVERTED_LINK_COLOR if row.S2 > row.E2 else DIRECT_LINK_COLOR
        circos.link((row.R, row.S1, row.E1), (row.Q, row.S2, row.E2), color=color, r1=r, r2=r)
    return len(links)
