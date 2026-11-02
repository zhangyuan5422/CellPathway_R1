# CellPathway R1

This repository contains the reproducible first stage of CellPathway:

1. enrichment of noncoding de novo mutations (DNMs) in two fixed enhancer universes;
2. TAD-constrained nomination of candidate genes for DNMs overlapping a selected enhancer set.

The packaged enhancer universes are fixed and validated:

- **EnhancerAtlas:** 52 callable/mappability-filtered tissue or cell-type enhancer sets;
- **scEnhancer:** 165 cell-type enhancer sets from the reproduced manuscript `part` universe (sets with source length under 100 Mb, excluding the added craniofacial set).

The code and packaged inputs were exported from the audited `reproduction_results` analysis. Coordinates use **GRCh38**.

## Repository layout

```text
data/
  dnvs/                  prepared noncoding DNM BED files at CADD >= 0, 10, 15, 20
  enhancers/
    atlas/               52 EnhancerAtlas BED files
    scenhancer/          165 scEnhancer BED files
    atlas_manifest.tsv
    scenhancer_manifest.tsv
  tad/                   brain and CNCC TADs plus the gene BigBed
scripts/
  validate_inputs.py
  run_enrichment.py
  annotate_tad.py
  run_step1.sh
results/                  generated locally; ignored by Git
```

BED inputs are individually gzip-compressed. The analysis scripts read `.bed.gz` directly.

## Installation

Python 3.10 or newer and `bedtools` are required.

```bash
conda create -n cellpathway-r1 -c conda-forge -c bioconda \
  python=3.10 bedtools numpy pandas scipy statsmodels pybigwig
conda activate cellpathway-r1
```

Alternatively, install `bedtools` separately and run `pip install -r requirements.txt`.

## Validate the packaged universes

```bash
python scripts/validate_inputs.py
```

Validation must report exactly 52 Atlas sets and 165 scEnhancer sets. It also checks file presence, interval counts, covered base pairs, and SHA-256 checksums against the manifests.

## Run DNM enrichment

Run both resources across all packaged cohorts and CADD thresholds:

```bash
bash scripts/run_step1.sh
```

Or run one resource:

```bash
python scripts/run_enrichment.py --resource atlas --workers 8
python scripts/run_enrichment.py --resource scenhancer --workers 8
```

Restrict a run to one threshold and cohort:

```bash
python scripts/run_enrichment.py \
  --resource atlas \
  --threshold 10 \
  --cohort spark \
  --workers 8
```

Results are written under `results/enrichment/<resource>/`. For each cohort and CADD threshold, the output reports enhancer length, DNM overlap count, DNM count in the remainder of the enhancer union, fold enrichment, one-sided Fisher exact P value, and Benjamini-Hochberg FDR across the complete resource universe.

This is the manuscript reproduction test: it uses enhancer base-pair length as the background and preserves the original interval-counting behavior. It should not be interpreted as a sequence-context- or mutation-rate-calibrated model.

## TAD-constrained candidate-gene annotation

After selecting a phenotype-relevant enriched enhancer set, annotate its overlapping DNMs. For example, fetal-brain enhancers and SPARK DNMs at CADD >= 10:

```bash
python scripts/annotate_tad.py \
  --dnv data/dnvs/cadd_10/dnm_spark.bed.gz \
  --resource atlas \
  --enhancer-set Fetal_brain \
  --tad brain \
  --output results/tad/spark_cadd10_fetal_brain.tsv
```

For craniofacial analyses, use `--tad cncc`. The output contains the overlapping enhancer interval, TAD interval, boundary-crossing status, and all protein-coding genes in that TAD.

TAD co-membership nominates genes in the same regulatory neighborhood; it does **not** establish a direct enhancer-target relationship.

## Input conventions

- Enhancer and TAD coordinates are GRCh38 BED coordinates.
- Packaged DNM BED files reproduce the original analysis and contain chromosome, position/start, and position/end without headers.
- Custom DNM inputs must be tab-delimited BED files with at least three columns and chromosome names matching the enhancer files (`chr1`, `chr2`, ...).
- To use custom DNMs, arrange them as `DIR/cadd_<threshold>/dnm_<cohort>.bed` or `.bed.gz` and pass `--dnv-root DIR`.

## Reproducibility provenance

- Atlas source: `reproduction_results/work/02_enhanceratlas/qc_mappable_hg38` (all 52 sets).
- scEnhancer source and selection: `reproduction_results/work/03_scenhancer/cluster_qc_hg38`, using the reproduced non-unique `part` rule and its asserted count of 165.
- Enrichment statistic: adapted from `reproduction_results/scripts/07_run_enrichment.py`.
- TAD resources: GSE77565 cortical-plate TADs for brain analyses and CNCC TADs for craniofacial analyses.

