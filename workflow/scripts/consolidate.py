#!/usr/bin/env python3
"""Consolidate NCBI metadata, BUSCO and antiSMASH results into per-assembly tables.

Every assembly in the NCBI metadata snapshot gets a row, whether or not its
BUSCO and antiSMASH runs have finished, so partial runs are visible as such.
Only small files are read (BUSCO summaries and antiSMASH region GenBank files),
so the script also works on a partial copy of the data folder.

Outputs (in --out-dir):
    BAGS_assemblies.tsv   one row per assembly
    BAGS_bgc_regions.tsv  one row per antiSMASH BGC region

Usage:
    python workflow/scripts/consolidate.py --data-dir data \
        --ncbi-report data/genomes/assembly_data_report.jsonl
"""

import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import Counter

# antiSMASH product classes that are characteristic of bacteria. In Rhizopus
# these point to endosymbiont (Mycetohabitans) or contaminant sequence.
BACTERIAL_TYPE_CLASSES = {
    "transAT-PKS", "transAT-PKS-like", "hserlactone", "butyrolactone",
    "phenazine", "lassopeptide", "ectoine", "thioamitides", "redox-cofactor",
    "phosphonate", "NRP-metallophore", "RiPP-like",
}
BACTERIAL_TYPE_PREFIXES = ("lanthipeptide",)

BUSCO_LINE = re.compile(
    r"C:(?P<C>[\d.]+)%\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],"
    r"F:(?P<F>[\d.]+)%,M:(?P<M>[\d.]+)%,n:(?P<n>\d+)"
)
ACCESSION = re.compile(r"^(GC[AF]_\d+\.\d+)")


def accession_of(name):
    match = ACCESSION.match(name)
    return match.group(1) if match else None


def snake_case(obj):
    """Recursively convert camelCase keys (datasets CLI jsonl) to snake_case (REST API)."""
    if isinstance(obj, dict):
        return {re.sub(r"(?<!^)(?=[A-Z])", "_", k).lower(): snake_case(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [snake_case(v) for v in obj]
    return obj


def read_ncbi_report(path):
    """Read the datasets CLI assembly_data_report.jsonl or a REST API dataset_report JSON."""
    with open(path) as handle:
        text = handle.read()
    try:
        data = json.loads(text)
        reports = data["reports"] if isinstance(data, dict) else data
    except json.JSONDecodeError:
        reports = [json.loads(line) for line in text.splitlines() if line.strip()]
    return [snake_case(r) for r in reports]


def ncbi_fields(report):
    organism = report.get("organism", {})
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    atypical = info.get("atypical", {})
    names = organism.get("infraspecific_names", {})
    organism_name = organism.get("organism_name", "")
    words = organism_name.split()
    variety = ""
    if "var." in words:
        variety = " ".join(words[words.index("var.") + 1:words.index("var.") + 2])
    return {
        "organism_name": organism_name,
        "species": " ".join(words[:2]),
        "variety": variety,
        "strain": names.get("strain") or names.get("isolate", ""),
        "assembly_level": info.get("assembly_level", ""),
        "is_atypical": bool(atypical.get("is_atypical", False)),
        "atypical_warnings": "; ".join(atypical.get("warnings", [])),
        "refseq_category": info.get("refseq_category", ""),
        "total_sequence_length": stats.get("total_sequence_length", ""),
        "contig_n50": stats.get("contig_n50", ""),
        "scaffold_n50": stats.get("scaffold_n50", ""),
        "number_of_contigs": stats.get("number_of_contigs", ""),
        "gc_percent": stats.get("gc_percent", ""),
        "release_date": info.get("release_date", ""),
        "submitter": info.get("submitter", ""),
    }


def read_busco(busco_dir):
    """Return BUSCO percentages from a BUSCO output folder, or None if not finished."""
    for path in sorted(glob.glob(os.path.join(busco_dir, "short_summary*.json"))):
        with open(path) as handle:
            data = json.load(handle)
        results = data["results"]
        return {
            "busco_C": results["Complete percentage"],
            "busco_S": results["Single copy percentage"],
            "busco_D": results["Multi copy percentage"],
            "busco_F": results["Fragmented percentage"],
            "busco_M": results["Missing percentage"],
            "busco_n": results["n_markers"],
            "busco_lineage": data.get("lineage_dataset", {}).get("name", ""),
        }
    for path in sorted(glob.glob(os.path.join(busco_dir, "short_summary*.txt"))):
        with open(path) as handle:
            match = BUSCO_LINE.search(handle.read())
        if match:
            values = {f"busco_{k}": float(v) for k, v in match.groupdict().items()}
            values["busco_n"] = int(values["busco_n"])
            values["busco_lineage"] = ""
            return values
    return None


def is_bacterial_type(product):
    return product in BACTERIAL_TYPE_CLASSES or product.startswith(BACTERIAL_TYPE_PREFIXES)


def read_region_gbk(path):
    """Parse one antiSMASH region GenBank file into a region record."""
    with open(path) as handle:
        text = handle.read()
    features, _, origin = text.partition("\nORIGIN")
    # The region feature runs until the next feature key at 5-space indent.
    block = re.search(r"\n     region .*?(?=\n     \S|\nORIGIN|$)", features, re.S)
    block = block.group(0) if block else ""
    products = sorted(set(re.findall(r'/product="([^"]+)"', block)))
    sequence = re.sub(r"[^acgtnACGTN]", "", origin).upper()
    acgt = sum(sequence.count(base) for base in "ACGT")
    start = re.search(r"Orig\. start\s*::\s*(\d+)", text)
    end = re.search(r"Orig\. end\s*::\s*(\d+)", text)
    number = re.search(r'/region_number="(\d+)"', block)
    edge = re.search(r'/contig_edge="(\w+)"', block)
    name = os.path.basename(path)
    return {
        "contig": name.rsplit(".region", 1)[0],
        "region_number": int(number.group(1)) if number else "",
        "orig_start": int(start.group(1)) if start else "",
        "orig_end": int(end.group(1)) if end else "",
        "length": len(sequence),
        "products": ";".join(products),
        "hybrid": len(products) > 1,
        "bacterial_type": any(is_bacterial_type(p) for p in products),
        "contig_edge": edge.group(1) == "True" if edge else "",
        "region_gc_percent": round(100 * (sequence.count("G") + sequence.count("C")) / acgt, 1) if acgt else "",
    }


def read_antismash(antismash_dir):
    """Return region records for a finished antiSMASH run, or None if not finished."""
    if not os.path.exists(os.path.join(antismash_dir, "index.html")):
        return None
    paths = sorted(glob.glob(os.path.join(antismash_dir, "*.region*.gbk")))
    return [read_region_gbk(p) for p in paths]


def find_output_dirs(parent):
    """Map accession -> output folder name for every result folder under parent."""
    if not os.path.isdir(parent):
        return {}
    return {accession_of(d): d for d in os.listdir(parent) if accession_of(d)}


def consolidate(data_dir, ncbi_report):
    reports = {r["accession"]: r for r in read_ncbi_report(ncbi_report)}
    busco_dirs = find_output_dirs(os.path.join(data_dir, "busco"))
    antismash_dirs = find_output_dirs(os.path.join(data_dir, "Antismash"))

    unknown = sorted((set(busco_dirs) | set(antismash_dirs)) - set(reports))
    if unknown:
        print(f"Warning: {len(unknown)} result folders are not in the NCBI report: "
              f"{', '.join(unknown)}", file=sys.stderr)

    assemblies, regions = [], []
    for accession in sorted(set(reports) | set(unknown)):
        report = reports.get(accession, {})
        # Name the assembly after its result folder, i.e. the .fna file stem.
        assembly = busco_dirs.get(accession) or antismash_dirs.get(accession)
        if not assembly:
            assembly_name = report.get("assembly_info", {}).get("assembly_name", "")
            assembly = f"{accession}_{assembly_name.replace(' ', '_')}_genomic"
        paired = report.get("paired_accession", "")
        row = {
            "assembly": assembly,
            "accession": accession,
            "in_ncbi_report": bool(report),
            # A GCF copy of a GCA assembly that is also in the set.
            "twin_of": paired if accession.startswith("GCF_") and paired in reports else "",
        }
        row.update(ncbi_fields(report) if report else {})

        busco = read_busco(os.path.join(data_dir, "busco", busco_dirs[accession])) if accession in busco_dirs else None
        row["busco_done"] = busco is not None
        row.update(busco or {})

        found = read_antismash(os.path.join(data_dir, "Antismash", antismash_dirs[accession])) if accession in antismash_dirs else None
        row["antismash_done"] = found is not None
        if found is not None:
            row["total_regions"] = len(found)
            row["hybrid_regions"] = sum(r["hybrid"] for r in found)
            row["bacterial_type_regions"] = sum(r["bacterial_type"] for r in found)
            classes = Counter(p for r in found for p in r["products"].split(";") if p)
            row.update({f"bgc_{cls}": n for cls, n in classes.items()})
            regions.extend({"assembly": assembly, "accession": accession, **r} for r in found)
        assemblies.append(row)
    return assemblies, regions


def write_tsv(path, rows, fixed_columns):
    extra = Counter(k for row in rows for k in row if k not in fixed_columns)
    columns = fixed_columns + [k for k, _ in extra.most_common()]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


ASSEMBLY_COLUMNS = [
    "assembly", "accession", "in_ncbi_report", "twin_of", "organism_name", "species", "variety",
    "strain", "assembly_level", "is_atypical", "atypical_warnings", "refseq_category",
    "total_sequence_length", "contig_n50", "scaffold_n50", "number_of_contigs", "gc_percent",
    "release_date", "submitter", "busco_done", "busco_C", "busco_S", "busco_D", "busco_F",
    "busco_M", "busco_n", "busco_lineage", "antismash_done", "total_regions", "hybrid_regions",
    "bacterial_type_regions",
]
REGION_COLUMNS = [
    "assembly", "accession", "contig", "region_number", "orig_start", "orig_end", "length",
    "products", "hybrid", "bacterial_type", "contig_edge", "region_gc_percent",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", default="data", help="BAGS data folder (default: data)")
    parser.add_argument("--ncbi-report", default="data/genomes/assembly_data_report.jsonl",
                        help="NCBI assembly report: datasets CLI .jsonl or REST API dataset_report JSON")
    parser.add_argument("--out-dir", default="data/consolidated", help="Output folder (default: data/consolidated)")
    args = parser.parse_args()

    assemblies, regions = consolidate(args.data_dir, args.ncbi_report)
    os.makedirs(args.out_dir, exist_ok=True)
    write_tsv(os.path.join(args.out_dir, "BAGS_assemblies.tsv"), assemblies, ASSEMBLY_COLUMNS)
    write_tsv(os.path.join(args.out_dir, "BAGS_bgc_regions.tsv"), regions, REGION_COLUMNS)

    print(f"{len(assemblies)} assemblies: "
          f"{sum(r['busco_done'] for r in assemblies)} with BUSCO, "
          f"{sum(r['antismash_done'] for r in assemblies)} with antiSMASH, "
          f"{sum(bool(r['twin_of']) for r in assemblies)} GCF twins; "
          f"{len(regions)} BGC regions -> {args.out_dir}")


if __name__ == "__main__":
    main()
