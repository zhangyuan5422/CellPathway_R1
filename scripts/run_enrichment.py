#!/usr/bin/env python3
"""Run the reproduced enhancer-length Fisher enrichment test."""

from __future__ import annotations

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

from validate_inputs import EXPECTED, load_manifest


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_THRESHOLDS = (0, 10, 15, 20)


def require_bedtools() -> str:
    executable = shutil.which("bedtools")
    if executable is None:
        raise RuntimeError("bedtools was not found on PATH")
    return executable


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def bed_length(path: Path) -> int:
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


def build_union(paths: list[Path], output: Path, bedtools: str) -> None:
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


def overlap_count(task: tuple[str, Path, Path]) -> int:
    bedtools, enhancer, dnv = task
    process = subprocess.Popen(
        [bedtools, "intersect", "-a", str(enhancer), "-b", str(dnv), "-c"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    total = 0
    assert process.stdout is not None
    for line in process.stdout:
        total += int(line.rstrip("\n").split("\t")[-1])
    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"bedtools failed for {enhancer} and {dnv}: {stderr}")
    return total


def discover_dnv(dnv_root: Path, threshold: int, cohort: str) -> Path:
    base = dnv_root / f"cadd_{threshold}" / f"dnm_{cohort}.bed"
    for candidate in (base, Path(f"{base}.gz")):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"No DNM BED found at {base} or {base}.gz")


def available_cohorts(dnv_root: Path, threshold: int) -> list[str]:
    folder = dnv_root / f"cadd_{threshold}"
    names = []
    for path in sorted(folder.glob("dnm_*.bed*")):
        name = path.name.removeprefix("dnm_").removesuffix(".gz").removesuffix(".bed")
        names.append(name)
    return list(dict.fromkeys(names))


def analyze_resource(
    resource: str,
    dnv_root: Path,
    output_root: Path,
    thresholds: list[int],
    requested_cohorts: list[str] | None,
    workers: int,
    check_inputs: bool,
) -> None:
    bedtools = require_bedtools()
    manifest = load_manifest(resource, check_content=check_inputs)
    if len(manifest) != EXPECTED[resource]:
        raise ValueError(f"{resource}: enhancer universe count changed")

    enhancer_root = ROOT / "data" / "enhancers" / resource
    names = manifest["enhancer_set"].tolist()
    paths = [enhancer_root / value for value in manifest["file"]]
    lengths = manifest["bp"].astype(int).tolist()
    resource_output = output_root / resource
    resource_output.mkdir(parents=True, exist_ok=True)
    union = resource_output / "enhancer_union.bed"
    build_union(paths, union, bedtools)
    union_length = bed_length(union)

    manifest[["enhancer_set", "intervals", "bp"]].to_csv(
        resource_output / "enhancer_lengths.tsv", sep="\t", index=False
    )
    summary_rows = []
    for threshold in thresholds:
        cohorts = requested_cohorts or available_cohorts(dnv_root, threshold)
        if not cohorts:
            raise ValueError(f"No cohorts found under {dnv_root / f'cadd_{threshold}'}")
        threshold_output = resource_output / f"cadd_{threshold}"
        threshold_output.mkdir(parents=True, exist_ok=True)
        for cohort in cohorts:
            dnv = discover_dnv(dnv_root, threshold, cohort)
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
            fdr = multipletests(p_values, method="fdr_bh")[1]
            frame = pd.DataFrame({
                "Enhancer_set": names,
                "Enhc_bp": lengths,
                f"dnm_overlap_{cohort}": overlaps,
                f"dnm_overlap_other_enhancers_{cohort}": [
                    union_overlap - value for value in overlaps
                ],
                f"Fold_enrichment_{cohort}": folds,
                f"P_value_{cohort}": p_values,
                f"P_FDR_{cohort}": fdr,
            })
            result_path = threshold_output / f"enhancer_enrichment_{cohort}_{threshold}.csv"
            frame.to_csv(result_path, index=False)
            summary_rows.append([
                resource, threshold, cohort, len(frame), union_length,
                union_overlap, result_path.relative_to(ROOT),
            ])
            print(
                f"{resource}: CADD>={threshold}, {cohort}: "
                f"{len(frame)} sets, {union_overlap} DNM overlaps in union"
            )

    with (resource_output / "analysis_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "resource", "cadd_threshold", "cohort", "enhancer_sets",
            "union_bp", "dnm_in_union", "result_file",
        ])
        writer.writerows(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--resource", choices=("atlas", "scenhancer", "all"), default="all"
    )
    parser.add_argument("--dnv-root", type=Path, default=ROOT / "data" / "dnvs")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "enrichment")
    parser.add_argument("--threshold", type=int, action="append", choices=DEFAULT_THRESHOLDS)
    parser.add_argument("--cohort", action="append")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--quick-input-check", action="store_true",
        help="Skip BED-content and checksum validation before analysis.",
    )
    args = parser.parse_args()
    resources = ("atlas", "scenhancer") if args.resource == "all" else (args.resource,)
    thresholds = args.threshold or list(DEFAULT_THRESHOLDS)
    for resource in resources:
        analyze_resource(
            resource=resource,
            dnv_root=args.dnv_root.resolve(),
            output_root=args.output.resolve(),
            thresholds=thresholds,
            requested_cohorts=args.cohort,
            workers=args.workers,
            check_inputs=not args.quick_input_check,
        )


if __name__ == "__main__":
    main()

