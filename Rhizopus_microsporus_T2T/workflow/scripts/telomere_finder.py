#!/usr/bin/env python3
"""
Find the telomeric repeat motif of an assembly with tidk.

  1. `tidk explore` lists candidate repeat motifs.
  2. The top N candidates by abundance are kept.
  3. `tidk search` counts each candidate in windows along every sequence.
  4. Each candidate is scored by how much of its signal sits in the first or
     last window of a sequence, times the fraction of sequences with such a
     terminal hit. A real telomeric repeat clusters at sequence ends; an
     abundant but unrelated repeat is spread across the whole sequence.
     (Abundance alone is not enough: in R. stolonifer PRFJ02 the telomeric
     motif AACCACAACCAC ranks only fourth by explore count.)
  5. `tidk plot` is run for the best candidate.

Writes <outdir>/<name>/<name>_candidate_scores.tsv with every candidate and
its score. The motif used for a circos telomere track is set in the config,
so this script is only needed to choose (or check) that motif.

Usage:
  telomere_finder.py --fasta genome.fasta --name Genome_name --outdir results/telomere_finder
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict


def explore(fasta, name, outdir, minimum, maximum, threshold, distance):
    out_tsv = os.path.join(outdir, f"{name}_exploration.tsv")
    cmd = [
        "tidk",
        "explore",
        "--minimum",
        str(minimum),
        "--maximum",
        str(maximum),
        "--threshold",
        str(threshold),
        "--distance",
        str(distance),
        fasta,
    ]
    print("+ " + " ".join(cmd) + f" > {out_tsv}", flush=True)
    with open(out_tsv, "w") as fh:
        subprocess.run(cmd, check=True, stdout=fh)
    return out_tsv


def top_motifs(exploration_tsv, top_n):
    with open(exploration_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        count_col = next(c for c in reader.fieldnames if c.startswith("count_repeat_runs"))
        rows = [
            (row["canonical_repeat_unit"], int(row[count_col]))
            for row in reader
            if row.get("canonical_repeat_unit")
        ]
    rows.sort(key=lambda r: r[1], reverse=True)
    return rows[:top_n]


def search(fasta, motif, name, outdir, window):
    label = f"{name}_{motif}"
    cmd = [
        "tidk",
        "search",
        "--string",
        motif,
        "--window",
        str(window),
        "--output",
        label,
        "--dir",
        outdir,
        fasta,
    ]
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    return os.path.join(outdir, f"{label}_telomeric_repeat_windows.tsv")


def score_windows_tsv(windows_tsv, terminal_hit_min_repeats):
    """Score one candidate motif from its tidk search windows.

    Returns (score, sequences with a terminal hit, sequences with any signal).
    score = fraction of all repeats that sit in the first or last window of a
    sequence x fraction of sequences with a terminal hit.
    """
    per_seq_rows = defaultdict(list)
    with open(windows_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            per_seq_rows[row["id"]].append(
                (
                    int(row["window"]),
                    int(row["forward_repeat_number"]),
                    int(row["reverse_repeat_number"]),
                )
            )

    total_signal = 0
    terminal_signal = 0
    n_seqs = 0
    n_seqs_terminal_hit = 0
    for rows in per_seq_rows.values():
        rows.sort(key=lambda r: r[0])
        signals = [fwd + rev for _, fwd, rev in rows]
        seq_total = sum(signals)
        if seq_total == 0:
            continue
        n_seqs += 1
        total_signal += seq_total
        terminal_signal += signals[0] + signals[-1]
        if signals[0] >= terminal_hit_min_repeats or signals[-1] >= terminal_hit_min_repeats:
            n_seqs_terminal_hit += 1

    if total_signal == 0 or n_seqs == 0:
        return 0.0, 0, 0
    terminal_fraction = terminal_signal / total_signal
    coverage = n_seqs_terminal_hit / n_seqs
    return terminal_fraction * coverage, n_seqs_terminal_hit, n_seqs


def plot(windows_tsv, out_prefix):
    cmd = ["tidk", "plot", "--tsv", windows_tsv, "-o", out_prefix]
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--fasta", required=True, help="Assembly FASTA")
    ap.add_argument(
        "--name", required=True, help="Genome name used for the output folder and file names"
    )
    ap.add_argument("--outdir", default="output", help="Results go in <outdir>/<name>/")
    ap.add_argument(
        "--top-n", type=int, default=10, help="Number of most abundant candidates to test"
    )
    ap.add_argument("--min", type=int, default=5, dest="minimum", help="tidk explore --minimum")
    ap.add_argument("--max", type=int, default=12, dest="maximum", help="tidk explore --maximum")
    ap.add_argument("--threshold", type=int, default=100, help="tidk explore --threshold")
    ap.add_argument("--distance", type=float, default=0.01, help="tidk explore --distance")
    ap.add_argument("--window", type=int, default=10000, help="tidk search --window")
    ap.add_argument(
        "--terminal-hit-min-repeats",
        type=int,
        default=5,
        help="Repeats needed in a first/last window to count as a terminal hit",
    )
    ap.add_argument(
        "--low-confidence-threshold",
        type=float,
        default=0.1,
        help="Warn when the best score is below this (no motif clusters at sequence ends)",
    )
    ap.add_argument(
        "--explore-only", action="store_true", help="Only run tidk explore and list the candidates"
    )
    args = ap.parse_args()

    genome_outdir = os.path.join(args.outdir, args.name)
    os.makedirs(genome_outdir, exist_ok=True)

    print(f"=== {args.name}: tidk explore ===")
    exploration_tsv = explore(
        args.fasta,
        args.name,
        genome_outdir,
        args.minimum,
        args.maximum,
        args.threshold,
        args.distance,
    )
    candidates = top_motifs(exploration_tsv, args.top_n)
    if not candidates:
        print(f"No candidate repeats found for {args.name} (try lowering --threshold).")
        sys.exit(1)

    print(f"Top {len(candidates)} candidates by abundance:")
    for motif, count in candidates:
        print(f"  {motif}\t{count}")

    if args.explore_only:
        return

    print(f"=== {args.name}: tidk search on {len(candidates)} candidates ===")
    results = []
    for motif, count in candidates:
        windows_tsv = search(args.fasta, motif, args.name, genome_outdir, args.window)
        score, n_hit, n_seqs = score_windows_tsv(windows_tsv, args.terminal_hit_min_repeats)
        results.append((motif, count, score, n_hit, n_seqs, windows_tsv))
        print(
            f"  {motif}\texplore_count={count}\tterminal_score={score:.3f}\t"
            f"contigs_with_terminal_hit={n_hit}/{n_seqs}"
        )

    results.sort(key=lambda r: r[2], reverse=True)
    best_motif, _, best_score, _, _, best_windows_tsv = results[0]

    summary_tsv = os.path.join(genome_outdir, f"{args.name}_candidate_scores.tsv")
    with open(summary_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "canonical_repeat_unit",
                "explore_count",
                "terminal_score",
                "contigs_with_terminal_hit",
                "n_contigs_with_signal",
                "is_best",
            ]
        )
        for motif, count, score, n_hit, n_seqs, _ in results:
            writer.writerow([motif, count, f"{score:.4f}", n_hit, n_seqs, motif == best_motif])

    if best_score < args.low_confidence_threshold:
        print(
            f"WARNING: best candidate for {args.name} only scored {best_score:.4f} "
            f"(< {args.low_confidence_threshold}) -- no candidate showed convincing "
            f"chromosome-end clustering, so '{best_motif}' is likely NOT the true "
            f"telomeric repeat (the assembly may simply not reach its telomeres). "
            f"Inspect {summary_tsv} and the plots manually, or rerun with a larger --top-n."
        )

    print(
        f"=== {args.name}: best candidate = {best_motif} (score={best_score:.3f}) -> tidk plot ==="
    )
    plot_prefix = os.path.join(genome_outdir, f"{args.name}_{best_motif}_plot")
    plot(best_windows_tsv, plot_prefix)

    print(f"Done. Summary: {summary_tsv}")
    print(f"Best motif windows TSV (for downstream circos telomere tracks): {best_windows_tsv}")


if __name__ == "__main__":
    main()
