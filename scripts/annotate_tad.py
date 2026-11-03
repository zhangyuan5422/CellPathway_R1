#!/usr/bin/env python3
"""Find genes in TADs containing DNM-overlapping enhancers."""

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pyBigWig

from validate_inputs import load_enhancers


ROOT = Path(__file__).resolve().parents[1]
TAD_FILES = {
    "brain": ROOT / "data" / "tad" / "GSE77565_CP_TAD_boundary.hg38.bed",
    "cncc": ROOT / "data" / "tad" / "CNCC_TADS_hg38_tad.bed",
}
GENES = ROOT / "data" / "tad" / "genes_w_noncoding.bb"


def enhancer_dnv_overlaps(enhancer, dnv):
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        raise RuntimeError("bedtools was not found on PATH")
    process = subprocess.Popen(
        [bedtools, "intersect", "-a", str(enhancer), "-b", str(dnv), "-wa", "-wb"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            raise ValueError("bedtools returned fewer than six columns")
        yield (
            fields[0], int(fields[1]), int(fields[2]),
            fields[-3], int(fields[-2]), int(fields[-1]),
        )
    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"bedtools intersect failed: {stderr}")


def interval_overlap(frame, chrom, start, end):
    # Reproduces bedtools' handling of the historical zero-width DNM records.
    query_start = start - 1 if start == end and start > 0 else start
    query_end = end + 1 if start == end else end
    return frame[
        (frame["chrom"] == chrom)
        & (frame["start"] < query_end)
        & (frame["end"] > query_start)
    ]


def locate_tad(chrom, start, end, domains, boundaries):
    overlapping_boundaries = interval_overlap(boundaries, chrom, start, end)
    if overlapping_boundaries.empty:
        overlapping_domains = interval_overlap(domains, chrom, start, end).copy()
        if overlapping_domains.empty:
            return None, "outside_TAD"
        overlapping_domains["overlap_bp"] = (
            overlapping_domains["end"].clip(upper=end)
            - overlapping_domains["start"].clip(lower=start)
        )
        selected = overlapping_domains.sort_values(
            ["overlap_bp", "start", "end"], ascending=[False, True, True]
        ).iloc[0]
        return (chrom, int(selected.start), int(selected.end)), "in_TAD"

    first_boundary = int(overlapping_boundaries["start"].min())
    last_boundary = int(overlapping_boundaries["end"].max())
    chrom_domains = domains[domains["chrom"] == chrom]
    upstream = chrom_domains[chrom_domains["end"] <= first_boundary].sort_values("end")
    downstream = chrom_domains[chrom_domains["start"] >= last_boundary].sort_values("start")
    if upstream.empty or downstream.empty:
        return None, "boundary_without_flanking_TADs"
    left = upstream.iloc[-1]
    right = downstream.iloc[0]
    return (chrom, int(left.start), int(right.end)), "cross_TAD_boundary"


def genes_in_region(bigbed, region):
    if region is None:
        return []
    try:
        entries = bigbed.entries(region[0], region[1], region[2]) or []
    except RuntimeError:
        return []
    genes = []
    for entry in entries:
        fields = entry[2].split("\t")
        if len(fields) > 1 and fields[1] == "1":
            genes.append(fields[0])
    return list(dict.fromkeys(genes))


def enhancer_path(resource, enhancer_set):
    enhancers = load_enhancers(resource, calculate_stats=False)
    match = enhancers[enhancers["enhancer_set"] == enhancer_set]
    if len(match) != 1:
        choices = ", ".join(enhancers["enhancer_set"].head(10))
        raise ValueError(
            f"Enhancer set {enhancer_set!r} is not in {resource}; first choices: {choices}"
        )
    return ROOT / "data" / "enhancers" / resource / match.iloc[0]["file"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dnv", type=Path, required=True)
    parser.add_argument("--resource", choices=("atlas", "scenhancer"), required=True)
    parser.add_argument("--enhancer-set", required=True)
    parser.add_argument("--tad", choices=tuple(TAD_FILES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    dnv_path = args.dnv.resolve()
    if not dnv_path.is_file():
        raise FileNotFoundError(dnv_path)
    selected_enhancer = enhancer_path(args.resource, args.enhancer_set)
    tad = pd.read_csv(
        TAD_FILES[args.tad], sep="\t", header=None,
        names=["chrom", "start", "end", "label"],
    )
    tad["start"] = tad["start"].astype(int)
    tad["end"] = tad["end"].astype(int)
    domains = tad[tad["label"].astype(str) == "0"][["chrom", "start", "end"]]
    boundaries = tad[tad["label"].astype(str) == "1"][["chrom", "start", "end"]]

    rows = []
    bigbed = pyBigWig.open(str(GENES))
    cache = {}
    for echrom, estart, eend, vchrom, vstart, vend in enhancer_dnv_overlaps(
        selected_enhancer, dnv_path
    ):
        enhancer_key = (echrom, estart, eend)
        if enhancer_key not in cache:
            region, status = locate_tad(
                echrom, estart, eend, domains, boundaries
            )
            cache[enhancer_key] = (region, status, genes_in_region(bigbed, region))
        region, status, genes = cache[enhancer_key]
        rows.append({
            "variant": f"{vchrom}:{vstart}-{vend}",
            "enhancer_resource": args.resource,
            "enhancer_set": args.enhancer_set,
            "overlapping_enhancer": f"{echrom}:{estart}-{eend}",
            "TAD_source": args.tad,
            "TAD": "" if region is None else f"{region[0]}:{region[1]}-{region[2]}",
            "TAD_status": status,
            "n_candidate_genes": len(genes),
            "candidate_genes": ",".join(genes),
        })
    bigbed.close()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "variant", "enhancer_resource", "enhancer_set", "overlapping_enhancer",
        "TAD_source", "TAD", "TAD_status", "n_candidate_genes", "candidate_genes",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(output, sep="\t", index=False)
    print(f"Wrote {len(rows)} enhancer-overlap/TAD annotations to {output}")


if __name__ == "__main__":
    main()
