#!/usr/bin/env python3
"""Draw the UFCG concatenated-gene tree with GSI values on the branches.

UFCG writes GSI (gene support index: how many of the single-copy core genes
support a split) as the internal node labels of the Newick file. The tree is
rooted on --outgroup, genome names are printed without underscores, and the
genomes listed in --highlight get an asterisk (*assembled in this study).
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--tree", required=True, help="UFCG Newick tree (e.g. concatenated_gsi_*.nwk)")
    ap.add_argument("--out", required=True, help="Output image")
    ap.add_argument("--outgroup", help="Leaf name to root on")
    ap.add_argument(
        "--highlight", nargs="*", default=[], help="Leaf names to mark with an asterisk"
    )
    ap.add_argument("--n-genes", type=int, help="Number of core genes (shown in the x-axis label)")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    tree = Phylo.read(args.tree, "newick")
    if args.outgroup:
        tree.root_with_outgroup(args.outgroup)
    tree.ladderize()

    def leaf_label(clade):
        if not clade.is_terminal():
            return None
        name = clade.name.replace("_", " ")
        return f"{name} *" if clade.name in args.highlight else name

    n_leaves = tree.count_terminals()
    fig, ax = plt.subplots(figsize=(10, 0.45 * n_leaves + 1.5))
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=leaf_label,
        label_colors=lambda name: "black",
    )

    # GSI next to each internal node (Phylo's own branch labels sit mid-branch
    # and overlap on the very short R. microsporus branches). Node positions
    # follow Phylo.draw: leaves at y = 1..n, internal nodes midway between
    # their first and last child, x = distance from the root.
    depths = tree.depths()
    y_pos = {leaf: i for i, leaf in enumerate(tree.get_terminals(), start=1)}
    for clade in tree.get_nonterminals(order="postorder"):
        y_pos[clade] = (y_pos[clade.clades[0]] + y_pos[clade.clades[-1]]) / 2
        if clade != tree.root and clade.confidence is not None:
            ax.text(
                depths[clade],
                y_pos[clade] - 0.08,
                f"{clade.confidence:g}",
                ha="right",
                va="bottom",
                fontsize=7,
                color="dimgrey",
            )
    for text in ax.texts:  # leaf names in italics
        if text.get_text().strip().startswith("Rhizopus"):
            text.set_fontstyle("italic")
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("left", "top", "right"):
        ax.spines[side].set_visible(False)
    genes = f"{args.n_genes} " if args.n_genes else ""
    ax.set_xlabel(f"Substitutions per site (concatenated {genes}UFCG core genes)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
