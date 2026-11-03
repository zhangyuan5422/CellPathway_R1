#!/usr/bin/env python3
"""Check the packaged enhancer BED files."""

import argparse
import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {"atlas": 52, "scenhancer": 165}


def bed_stats(path):
    opener = gzip.open if path.suffix == ".gz" else open
    rows = 0
    bp = 0
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: fewer than three BED columns")
            start, end = int(fields[1]), int(fields[2])
            if end < start:
                raise ValueError(f"{path}:{line_number}: end is smaller than start")
            rows += 1
            bp += end - start
    return rows, bp


def enhancer_name(path):
    name = path.name.removesuffix(".gz").removesuffix(".bed")
    return name.removesuffix(".hg38")


def load_enhancers(resource, calculate_stats=True):
    if resource not in EXPECTED:
        raise ValueError(f"Unknown enhancer resource: {resource}")
    base = ROOT / "data" / "enhancers" / resource
    paths = sorted(
        set(base.glob("*.bed")) | set(base.glob("*.bed.gz")),
        key=lambda path: path.name,
    )
    rows = []
    for set_index, path in enumerate(paths, 1):
        intervals, bp = bed_stats(path) if calculate_stats else (None, None)
        rows.append({
            "set_index": set_index,
            "enhancer_set": enhancer_name(path),
            "file": path.name,
            "intervals": intervals,
            "bp": bp,
        })
    frame = pd.DataFrame(rows)
    if len(frame) != EXPECTED[resource]:
        raise ValueError(
            f"{resource}: expected {EXPECTED[resource]} sets, found {len(frame)}"
        )
    if frame["enhancer_set"].duplicated().any() or frame["file"].duplicated().any():
        raise ValueError(f"{resource}: duplicate enhancer set or file")
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quick", action="store_true",
        help="Check enhancer directory file counts without parsing BED contents.",
    )
    args = parser.parse_args()
    for resource in EXPECTED:
        frame = load_enhancers(resource, calculate_stats=not args.quick)
        message = f"{resource}: PASS; {len(frame)} sets"
        if not args.quick:
            message += (
                f"; {int(frame['intervals'].sum()):,} intervals; "
                f"{int(frame['bp'].sum()):,} summed bp"
            )
        print(message)


if __name__ == "__main__":
    main()
