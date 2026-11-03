#!/usr/bin/env python3
"""Run the enhancer-length Fisher enrichment test."""

import argparse
import csv
import gzip
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from validate_inputs import EXPECTED, load_enhancers


ROOT = Path(__file__).resolve().parents[1]


def require_bedtools():
    bedtools = shutil.which("bedtools")
    if not bedtools:
        raise RuntimeError("bedtools was not found on PATH")
    return bedtools


def open_text(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def bed_length(path):
    total = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: fewer than three BED columns")
            start, end = int(fields[1]), int(fields[2])
            if end < start:
                raise ValueError(f"{path}:{line_number}: end is smaller than start")
            total += end - start
    return total


def build_union(paths, output, bedtools):
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cellpathway_union_") as temp_name:
        temp = Path(temp_name)
        combined = temp / "combined.bed"
        sorted_bed = temp / "sorted.bed"
        with combined.open("w") as destination:
            for path in paths:
                with open_text(path) as source:
                    shutil.copyfileobj(source, destination)
        with sorted_bed.open("w") as destination:
            subprocess.run(
                ["sort", "-k1,1", "-k2,2n", str(combined)],
                check=True, stdout=destination,
            )
        with output.open("w") as destination:
            subprocess.run(
                [bedtools, "merge", "-i", str(sorted_bed)],
                check=True, stdout=destination,
            )


def overlap_count(task):
    bedtools, enhancer, dnv = task
    process = subprocess.Popen(
        [bedtools, "intersect", "-a", str(enhancer), "-b", str(dnv), "-c"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    total = 0
    for line in process.stdout:
        total += int(line.rstrip("\n").split("\t")[-1])
    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"bedtools failed for {enhancer} and {dnv}: {stderr}")
    return total


def analyze_resource(resource, dnv, output_root, workers):
    bedtools = require_bedtools()
    enhancers = load_enhancers(resource, calculate_stats=True)
    if len(enhancers) != EXPECTED[resource]:
        raise ValueError(f"{resource}: enhancer universe count changed")

    enhancer_root = ROOT / "data" / "enhancers" / resource
    names = enhancers["enhancer_set"].tolist()
    paths = [enhancer_root / value for value in enhancers["file"]]
    lengths = enhancers["bp"].astype(int).tolist()
    resource_output = output_root / resource
    resource_output.mkdir(parents=True, exist_ok=True)
    union = resource_output / "enhancer_union.bed"
    build_union(paths, union, bedtools)
    union_length = bed_length(union)

    enhancers[["enhancer_set", "intervals", "bp"]].to_csv(
        resource_output / "enhancer_lengths.tsv", sep="\t", index=False
    )
    tasks = [(bedtools, path, dnv) for path in paths]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        overlaps = list(pool.map(overlap_count, tasks))
    union_overlap = overlap_count((bedtools, union, dnv))

    folds = []
    p_values = []
    for overlap, enhancer_bp in zip(overlaps, lengths):
        other_dnv = union_overlap - overlap
        other_bp = union_length - enhancer_bp
        if enhancer_bp <= 0 or other_bp <= 0:
            raise ValueError(f"Invalid opportunity counts for {resource}")
        if other_dnv == 0:
            fold = np.inf if overlap else np.nan
        else:
            fold = (overlap / enhancer_bp) / (other_dnv / other_bp)
        table = [
            [int(overlap), int(enhancer_bp - overlap)],
            [int(other_dnv), int(other_bp - other_dnv)],
        ]
        _, p_value = fisher_exact(table, alternative="greater")
        folds.append(fold)
        p_values.append(p_value)

    frame = pd.DataFrame({
        "Enhancer_set": names,
        "Enhc_bp": lengths,
        "dnm_overlap": overlaps,
        "dnm_overlap_other_enhancers": [
            union_overlap - value for value in overlaps
        ],
        "Fold_enrichment": folds,
        "P_value": p_values,
        "P_FDR": multipletests(p_values, method="fdr_bh")[1],
    })
    result_path = resource_output / "enhancer_enrichment.csv"
    frame.to_csv(result_path, index=False)
    print(
        f"{resource}: {len(frame)} sets, "
        f"{union_overlap} DNM overlaps in union"
    )

    with (resource_output / "analysis_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "resource", "enhancer_sets", "union_bp",
            "dnm_in_union", "result_file",
        ])
        writer.writerow([
            resource, len(frame), union_length, union_overlap,
            result_path.relative_to(ROOT),
        ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--resource", choices=("atlas", "scenhancer", "all"), default="all"
    )
    parser.add_argument(
        "--dnv", type=Path, default=ROOT / "example" / "example.bed"
    )
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "enrichment")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    resources = ("atlas", "scenhancer") if args.resource == "all" else (args.resource,)
    dnv = args.dnv.resolve()
    if not dnv.is_file():
        parser.error(f"DNM BED file not found: {dnv}")
    for resource in resources:
        analyze_resource(
            resource=resource,
            dnv=dnv,
            output_root=args.output.resolve(),
            workers=args.workers,
        )


if __name__ == "__main__":
    main()
