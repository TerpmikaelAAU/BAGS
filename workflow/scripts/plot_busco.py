#!/usr/bin/env python3
"""Violin plot of BUSCO categories per assembly, with optional filters and grouping.

Reads BAGS_assemblies.tsv from consolidate.py. Only assemblies with a finished
BUSCO run are plotted, and GCF twins are always left out. Every violin gets one
jittered dot per assembly and a median line.

Examples:
    python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv
    python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv \
        --group-by species --exclude-atypical-warning contaminated --out data/consolidated/busco_species
    python workflow/scripts/plot_busco.py data/consolidated/BAGS_assemblies.tsv \
        --compare-unfiltered --exclude-atypical-warning contaminated --min-busco-c 90
"""

import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CATEGORIES = [
    ("busco_C", "Complete"),
    ("busco_S", "Single-copy"),
    ("busco_D", "Duplicated"),
    ("busco_F", "Fragmented"),
    ("busco_M", "Missing"),
]
# Categorical slots in fixed order; a group keeps its colour whatever is filtered out.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER_COLOR = "#898781"
INK, INK_SECONDARY, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"


def read_table(path):
    table = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])
    table["busco_done"] = table["busco_done"].astype(str) == "True"
    return table[table["twin_of"].isna()]


def apply_filters(table, args):
    """Return the filtered table and a human-readable description of each filter."""
    kept, notes = table, []
    for warning in args.exclude_atypical_warning or []:
        flagged = kept["atypical_warnings"].fillna("").str.contains(warning, case=False, regex=False)
        kept = kept[~flagged]
        notes.append(f"NCBI atypical '{warning}' excluded")
    if args.min_busco_c is not None:
        kept = kept[kept["busco_C"] >= args.min_busco_c]
        notes.append(f"BUSCO C ≥ {args.min_busco_c:g}%")
    if args.species:
        kept = kept[kept["species"].isin(args.species)]
        notes.append("species: " + ", ".join(args.species))
    return kept, notes


def assign_groups(table, full, column, min_size):
    """Label each assembly with its group; groups smaller than min_size become 'Other'.

    Group order and colours come from the full (unfiltered) table, so filtering
    never repaints the groups that remain.
    """
    counts = full[column].fillna("Unknown").value_counts()
    large = [g for g in counts.index if counts[g] >= min_size]
    colors = {g: PALETTE[i] for i, g in enumerate(large[:len(PALETTE)])}
    labels = table[column].fillna("Unknown").where(lambda s: s.isin(colors), "Other")
    order = [g for g in large if g in colors] + ["Other"]
    colors["Other"] = OTHER_COLOR
    return labels, order, colors


def draw(groups, order, colors, title, subtitle, italic_labels, out):
    order = [g for g in order if g in groups and len(groups[g])]
    n_groups = len(order)
    slot = 0.8 / n_groups
    fig, ax = plt.subplots(figsize=(1.5 + len(CATEGORIES) * (0.7 + 0.35 * n_groups), 4.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    rng = np.random.default_rng(0)

    for gi, group in enumerate(order):
        data = groups[group]
        color = colors[group]
        for ci, (column, _) in enumerate(CATEGORIES):
            values = data[column].to_numpy(dtype=float)
            x = ci - 0.4 + slot * (gi + 0.5)
            if len(values) >= 2 and np.ptp(values) > 0:
                parts = ax.violinplot(values, positions=[x], widths=slot * 0.9,
                                      showextrema=False, bw_method=0.3)
                for body in parts["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor(color)
                    body.set_alpha(0.28)
                    body.set_linewidth(1)
            jitter = rng.uniform(-slot * 0.22, slot * 0.22, len(values))
            ax.scatter(x + jitter, values, s=9, color=color, alpha=0.85, linewidths=0, zorder=3)
            median = np.median(values)
            ax.hlines(median, x - slot * 0.3, x + slot * 0.3, color=INK, linewidth=1.5, zorder=4)

    ax.set_xticks(range(len(CATEGORIES)))
    ax.set_xticklabels([label for _, label in CATEGORIES], color=INK)
    ax.set_xlim(-0.5, len(CATEGORIES) - 0.5)
    ax.set_ylim(-2, 102)
    ax.set_ylabel("BUSCOs (%)", color=INK_SECONDARY)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(colors=INK_SECONDARY, length=0)

    ax.set_title(title, loc="left", color=INK, fontsize=12, pad=22)
    ax.text(0, 1.02, subtitle, transform=ax.transAxes, color=INK_SECONDARY, fontsize=9)

    if n_groups > 1:
        handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=7, color=colors[g],
                              label=f"{g} (n={len(groups[g])})") for g in order]
        legend = ax.legend(handles=handles, frameon=False, loc="upper left",
                           bbox_to_anchor=(1.0, 1.0), labelcolor=INK)
        if italic_labels:
            for text, group in zip(legend.get_texts(), order):
                if group not in ("Other", "Unknown"):
                    text.set_fontstyle("italic")

    fig.tight_layout()
    for extension in ("png", "pdf"):
        fig.savefig(f"{out}.{extension}", dpi=300, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table", help="BAGS_assemblies.tsv from consolidate.py")
    parser.add_argument("--out", default="busco_violin", help="Output path without extension (default: busco_violin)")
    parser.add_argument("--group-by", help="Column to split violins by, e.g. species or assembly_level")
    parser.add_argument("--min-group-size", type=int, default=5,
                        help="Groups smaller than this are pooled as 'Other' (default: 5)")
    parser.add_argument("--compare-unfiltered", action="store_true",
                        help="Show all assemblies next to the filtered set instead of grouping")
    filters = parser.add_argument_group("filters")
    filters.add_argument("--exclude-atypical-warning", action="append", metavar="TEXT",
                         help="Drop assemblies whose NCBI atypical warning contains TEXT, "
                              "e.g. 'contaminated' (repeatable)")
    filters.add_argument("--min-busco-c", type=float, help="Drop assemblies with BUSCO Complete below this %%")
    filters.add_argument("--species", action="append", help="Keep only this species (repeatable)")
    args = parser.parse_args()

    if args.group_by and args.compare_unfiltered:
        parser.error("use either --group-by or --compare-unfiltered, not both")

    full = read_table(args.table)
    pending = int((~full["busco_done"]).sum())
    full = full[full["busco_done"]]
    filtered, notes = apply_filters(full, args)
    if filtered.empty:
        sys.exit("No assemblies left after filtering.")

    if args.compare_unfiltered:
        groups = {"All assemblies": full, "Filtered": filtered}
        order = list(groups)
        colors = {"All assemblies": OTHER_COLOR, "Filtered": PALETTE[0]}
        plotted = pd.concat([full.assign(group="All assemblies"), filtered.assign(group="Filtered")])
    elif args.group_by:
        labels, order, colors = assign_groups(filtered, full, args.group_by, args.min_group_size)
        groups = {g: filtered[labels == g] for g in order}
        plotted = filtered.assign(group=labels)
    else:
        groups, order, colors = {"All": filtered}, ["All"], {"All": PALETTE[0]}
        plotted = filtered.assign(group="All")

    subtitle = f"{len(filtered)} assemblies"
    if notes:
        subtitle += " · " + " · ".join(notes)
    if pending:
        subtitle += f" · {pending} without BUSCO results not shown"
    lineage = full["busco_lineage"].dropna().unique()
    title = "BUSCO completeness" + (f" ({', '.join(lineage)})" if len(lineage) else "")

    draw(groups, order, colors, title, subtitle, args.group_by == "species", args.out)

    columns = ["assembly", "accession", "species", "group"] + [c for c, _ in CATEGORIES] + ["busco_n"]
    plotted[columns].to_csv(f"{args.out}.tsv", sep="\t", index=False)
    print(f"Plotted {len(filtered)} assemblies -> {args.out}.png/.pdf/.tsv")


if __name__ == "__main__":
    main()
