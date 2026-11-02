#!/usr/bin/env python3
"""Validate the fixed CellPathway enhancer universes against their manifests."""

from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {"atlas": 52, "scenhancer": 165}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bed_stats(path: Path) -> tuple[int, int]:
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


def load_manifest(resource: str, check_content: bool = True) -> pd.DataFrame:
    manifest_path = ROOT / "data" / "enhancers" / f"{resource}_manifest.tsv"
    frame = pd.read_csv(manifest_path, sep="\t", dtype={"sha256": str})
    required = {"set_index", "enhancer_set", "file", "intervals", "bp", "sha256"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{manifest_path}: missing columns {sorted(missing)}")
    if len(frame) != EXPECTED[resource]:
        raise ValueError(
            f"{resource}: expected {EXPECTED[resource]} sets, found {len(frame)}"
        )
    if frame["enhancer_set"].duplicated().any() or frame["file"].duplicated().any():
        raise ValueError(f"{resource}: duplicate enhancer set or file in manifest")
    if frame["set_index"].tolist() != list(range(1, len(frame) + 1)):
        raise ValueError(f"{resource}: set_index must be consecutive and one-based")

    base = ROOT / "data" / "enhancers" / resource
    for row in frame.itertuples(index=False):
        path = base / row.file
        if not path.is_file():
            raise FileNotFoundError(path)
        if check_content:
            intervals, bp = bed_stats(path)
            if intervals != int(row.intervals) or bp != int(row.bp):
                raise ValueError(
                    f"{path}: manifest says {row.intervals} rows/{row.bp} bp, "
                    f"observed {intervals} rows/{bp} bp"
                )
            observed_hash = sha256(path)
            if observed_hash != row.sha256:
                raise ValueError(f"{path}: SHA-256 mismatch")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quick", action="store_true",
        help="Check manifests and file presence without rereading all BED contents.",
    )
    args = parser.parse_args()
    for resource in EXPECTED:
        frame = load_manifest(resource, check_content=not args.quick)
        print(
            f"{resource}: PASS; {len(frame)} sets; "
            f"{int(frame['intervals'].sum()):,} intervals; {int(frame['bp'].sum()):,} summed bp"
        )


if __name__ == "__main__":
    main()

