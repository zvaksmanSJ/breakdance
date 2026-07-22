# Breakdance

**Breakdance** is a starter ONT structural variant interpretation pipeline for long-read structural variant callsets.

It is designed to work from structural variant VCF outputs produced by:

- **Sniffles2**
- **cuteSV**
- **SVIM**

and convert them into a first-pass interpretation framework with:

- harmonized breakpoint junctions
- junction clustering
- optional interval annotation
- interpreted candidate events
- event prioritization
- lightweight breakpoint graph export
- cohort summary outputs

---

# What Breakdance does

For each sample, Breakdance:

1. reads SV caller outputs from Sniffles2, cuteSV, and/or SVIM
2. parses raw VCF records into a shared schema
3. adds local breakpoint-density hints
4. harmonizes similar calls across callers into consensus junctions
5. clusters nearby junctions into candidate event groups
6. optionally annotates breakpoints with:
   - genes
   - exons
   - promoters
   - enhancers
7. creates interpreted event candidates such as:
   - fusion candidates
   - exon disruptions
   - promoter hijack candidates
   - enhancer hijack candidates
   - interchromosomal rearrangements
   - complex rearrangement clusters
   - chromothripsis-like clusters
8. assigns simple event priority scores and tiers
9. exports a lightweight breakpoint graph
10. writes per-sample and cohort-level summary files

---

# Current scope

Breakdance is currently a **Version 1 interpretation framework**.

It is intended for:

- ONT SV review
- candidate prioritization
- long-read rearrangement triage
- multi-caller harmonization
- exploratory cancer SV interpretation

It is **not yet**:

- a full copy-number-aware reconstruction engine
- a JaBbA replacement
- a clinically validated interpretation system
- a formal ecDNA / BFB classifier
- a full chromothripsis inference engine

---

# Repository structure

Core modules:

- `schema.py`  
  Core dataclasses used across the pipeline.

- `vcf_parser.py`  
  Reads and normalizes VCFs from Sniffles2, cuteSV, and SVIM.

- `harmonize.py`  
  Harmonizes raw calls across callers into consensus breakpoint junctions.

- `clustering.py`  
  Clusters harmonized junctions into candidate event groups.

- `annotate.py`  
  Adds BED4-based interval annotations.

- `interpret.py`  
  Generates first-pass interpreted event candidates.

- `score.py`  
  Assigns event priority scores and tiers.

- `graph.py`  
  Builds a lightweight breakpoint graph object.

- `report.py`  
  Writes TSV and JSON outputs.

- `run_harmonize.py`  
  Main CLI entrypoint for the pipeline.

- `setup_references.py`  
  Creates a build-specific reference directory and generates annotation BEDs.

Optional / future:

- `dashboard_app.py`
- copy-number integration modules
- complex-event formalization modules
- graph analysis extensions

---

# Prerequisites

## 1. Python

Recommended:

- **Python 3.10+**
- Python 3.11 is a good choice

Example conda environment:

```bash
conda create -n ont_sv python=3.11 -y
conda activate ont_sv
```

---

## 2. Required Breakdance files

To run `run_harmonize.py`, you need these files present in the same directory:

- `run_harmonize.py`
- `schema.py`
- `vcf_parser.py`
- `harmonize.py`
- `clustering.py`
- `annotate.py`
- `interpret.py`
- `score.py`
- `graph.py`
- `report.py`

If you are using reference setup as well:

- `setup_references.py`

---

## 3. Python package requirements

### Required for the core pipeline
The current core Breakdance pipeline uses only Python standard library modules.

So for `run_harmonize.py`, you do **not** need to install extra Python packages beyond Python itself.

Used standard modules include:
- `argparse`
- `csv`
- `gzip`
- `json`
- `pathlib`
- `re`
- `shutil`
- `statistics`
- `sys`
- `urllib.request`

### Optional packages
Only needed if you later use the dashboard:

```bash
pip install streamlit pandas
```

---

## 4. Input BAM list

Breakdance expects a BAM list file with **one BAM path per line**.

Example:

```text
/path/to/sampleA.bam
/path/to/sampleB.bam
/path/to/sampleC.bam
```

Rules:

- blank lines are ignored
- lines starting with `#` are ignored
- all listed BAMs must exist

Breakdance uses each BAM path to determine:
- sample name
- where the corresponding SV caller outputs should be located

---

## 5. Existing SV caller outputs

`run_harmonize.py` does **not** run Sniffles2, cuteSV, or SVIM for you.

It expects their VCF outputs to already exist.

For each sample, it looks for:

### Sniffles2
- `sniffles2/sniffles2.vcf`
- `sniffles2/sniffles2.vcf.gz`

### cuteSV
- `cutesv/cutesv.vcf`
- `cutesv/cutesv.vcf.gz`

### SVIM
- `svim/variants.vcf`
- `svim/variants.vcf.gz`

---

## 6. Expected sample output directory layout

Breakdance supports two ways to find caller outputs.

### Option A: use `--output-root`

If you run with:

```bash
--output-root /path/to/all_sv_calls
```

and the sample name is `sample1`, Breakdance expects:

```text
/path/to/all_sv_calls/sample1/
  sniffles2/sniffles2.vcf(.gz)
  cutesv/cutesv.vcf(.gz)
  svim/variants.vcf(.gz)
```

### Option B: no `--output-root`

If `--output-root` is omitted, Breakdance expects outputs next to each BAM:

```text
<bam_parent>/<sample>_sv_callers/
```

Example:

BAM:
```text
/data/run1/sample1.bam
```

Expected caller output root:
```text
/data/run1/sample1_sv_callers/
```

---

## 7. Optional annotation files

Annotation files are optional, but strongly recommended if you want gene-level interpretation.

Supported flags:

- `--genes-bed`
- `--exons-bed`
- `--promoters-bed`
- `--enhancers-bed`

These should be **BED4** files:

```text
chrom    start    end    name
```

Examples:
```text
chr8    127735433    127742951    MYC
chr14   105586337    105864198    IGH
```

---

# Setting up references with `setup_references.py`

## Purpose

`setup_references.py` creates a standardized Breakdance reference directory for:

- `hg19`
- `hg38`
- `t2t`

It installs and organizes:

- genome FASTA
- annotation GTF
- BED4 files for:
  - all genes
  - all exons
  - all promoters
- optional enhancer BED
- a manifest file

This is the recommended way to prepare reference annotations for Breakdance.

---

## Supported builds

- `hg19`
- `hg38`
- `t2t`

Build aliases accepted by the script include:

- `hg19`, `grch37`, `37`
- `hg38`, `grch38`, `38`
- `t2t`, `hs1`, `chm13`

---

## Reference directory layout

If you run:

```bash
python setup_references.py --build hg38 --reference-root references
```

you will get:

```text
references/
  hg38/
    genome/
      genome.fa.gz
    annotation/
      gtf/
        genes.gtf.gz
      bed/
        genes.bed
        exons.bed
        promoters.bed
        enhancers.bed        # only if supplied
    source_files/
    manifests/
      reference_manifest.tsv
```

---

## What `setup_references.py` generates

From the full annotation GTF, it builds:

### 1. `genes.bed`
All genes in BED4 format.

Name field format:
```text
gene_name|gene_id|gene_type
```

### 2. `exons.bed`
All exons in BED4 format.

Name field format:
```text
gene_name|gene_id|transcript_id|exon_number
```

### 3. `promoters.bed`
Promoter intervals generated from transcript TSS positions.

Default promoter window:
- **2000 bp upstream**
- **500 bp downstream**

Name field format:
```text
gene_name|gene_id|transcript_id|strand
```

### 4. `enhancers.bed`
Optional only. This must be supplied separately because enhancer sources are not universal.

---

## Default behavior by build

### hg38
Uses default Gencode GRCh38 URLs for:
- FASTA
- GTF

### hg19
Uses default Gencode lift37 / GRCh37-mapped resources for:
- FASTA
- GTF

Operationally this is treated as `hg19`.

### t2t
No hard-coded default URLs are enforced for T2T because source conventions vary more.

For T2T, you should usually provide:
- `--fasta` or `--fasta-url`
- `--gtf` or `--gtf-url`

---

## Reference setup examples

### hg38 using defaults

```bash
python setup_references.py \
  --build hg38 \
  --reference-root references
```

### hg19 using defaults

```bash
python setup_references.py \
  --build hg19 \
  --reference-root references
```

### T2T using explicit URLs

```bash
python setup_references.py \
  --build t2t \
  --reference-root references \
  --fasta-url https://your-source.example/hs1.fa.gz \
  --gtf-url https://your-source.example/hs1.annotation.gtf.gz
```

### T2T using local files

```bash
python setup_references.py \
  --build t2t \
  --reference-root references \
  --fasta /path/to/hs1.fa.gz \
  --gtf /path/to/hs1.annotation.gtf.gz
```

### Using local hg38 files

```bash
python setup_references.py \
  --build hg38 \
  --reference-root references \
  --fasta /path/to/hg38.fa.gz \
  --gtf /path/to/gencode.v47.annotation.gtf.gz
```

### Install an enhancer BED

```bash
python setup_references.py \
  --build hg38 \
  --reference-root references \
  --enhancers-bed /path/to/enhancers.bed
```

---

## Optional `setup_references.py` arguments

### Required

- `--build`
- `--reference-root`

### Optional input sources

- `--fasta`
- `--gtf`
- `--enhancers-bed`

- `--fasta-url`
- `--gtf-url`
- `--enhancers-bed-url`

### Optional behavior flags

- `--symlink-local-files`  
  Symlink local files instead of copying them.

- `--promoter-upstream-bp`  
  Default: `2000`

- `--promoter-downstream-bp`  
  Default: `500`

- `--force`  
  Allow overwriting existing generated files.

---

## Using generated references in Breakdance

Example using hg38 references created by `setup_references.py`:

```bash
python run_harmonize.py \
  --bam-list bams_to_run.txt \
  --genes-bed references/hg38/annotation/bed/genes.bed \
  --exons-bed references/hg38/annotation/bed/exons.bed \
  --promoters-bed references/hg38/annotation/bed/promoters.bed \
  --enhancers-bed references/hg38/annotation/bed/enhancers.bed \
  --summary-dir harmonized_summary
```

If no enhancer BED was installed, omit `--enhancers-bed`.

---

# Running Breakdance

## Minimal run

```bash
python run_harmonize.py \
  --bam-list /path/to/bams_to_run.txt \
  --summary-dir harmonized_summary
```

## Run with shared caller output root

```bash
python run_harmonize.py \
  --bam-list /path/to/bams_to_run.txt \
  --output-root /path/to/all_sv_calls \
  --summary-dir harmonized_summary
```

## Run with annotation references

```bash
python run_harmonize.py \
  --bam-list /path/to/bams_to_run.txt \
  --output-root /path/to/all_sv_calls \
  --genes-bed references/hg38/annotation/bed/genes.bed \
  --exons-bed references/hg38/annotation/bed/exons.bed \
  --promoters-bed references/hg38/annotation/bed/promoters.bed \
  --summary-dir harmonized_summary
```

## Run with focus genes

```bash
python run_harmonize.py \
  --bam-list /path/to/bams_to_run.txt \
  --output-root /path/to/all_sv_calls \
  --genes-bed references/hg38/annotation/bed/genes.bed \
  --exons-bed references/hg38/annotation/bed/exons.bed \
  --promoters-bed references/hg38/annotation/bed/promoters.bed \
  --focus-genes MYC,BCL2,TP53 \
  --summary-dir harmonized_summary
```

---

# `run_harmonize.py` arguments

## Required

- `--bam-list`  
  Text file with one BAM path per line.

## Optional

- `--output-root`  
  Common root directory containing per-sample caller outputs.

- `--tolerance-bp`  
  Breakpoint tolerance for harmonization.  
  Default: `100`

- `--cluster-distance-bp`  
  Breakpoint distance threshold for clustering harmonized junctions.  
  Default: `10000`

- `--summary-dir`  
  Directory for the cohort summary TSV.  
  Default: current directory

- `--genes-bed`
- `--exons-bed`
- `--promoters-bed`
- `--enhancers-bed`

- `--focus-genes`  
  Comma-separated gene list to boost event ranking.

---

# Per-sample outputs

For each sample, Breakdance writes:

## 1. Harmonized breakpoint TSV

```text
<sample>.harmonized_breakpoints.tsv
```

Contains:
- consensus junctions
- breakpoint coordinates
- SV type
- supporting callers
- support summaries
- confidence
- annotations
- source VCF paths

## 2. Junction clusters TSV

```text
<sample>.junction_clusters.tsv
```

Contains:
- cluster IDs
- member junction IDs
- chromosomes involved
- interchromosomal content
- cluster notes
- aggregated annotations

## 3. Interpreted events TSV

```text
<sample>.interpreted_events.tsv
```

Contains:
- interpreted event IDs
- event type
- linked junctions / cluster IDs
- genes involved
- confidence
- priority score
- priority tier
- event details

## 4. Breakpoint graph JSON

```text
<sample>.breakpoint_graph.json
```

Contains:
- graph nodes
- graph edges
- metadata

## 5. Dashboard JSON

```text
<sample>.dashboard.json
```

Contains a bundled JSON summary of:
- sample-level counts
- junctions
- clusters
- events
- graph

Even if you do not use a dashboard, this file can still be useful for downstream custom visualization.

---

# Cohort output

## Cohort summary TSV

```text
harmonized_cohort_summary.tsv
```

Contains one row per sample with:

- sample name
- BAM path
- sample output directory
- callers found
- raw call count
- harmonized junction count
- cluster count
- event count
- largest cluster size
- paths to sample-level outputs

---

# Interpretation logic

## Junction-level event types

Breakdance currently emits:

- `fusion_candidate`
- `exon_disruption`
- `promoter_hijack_candidate`
- `enhancer_hijack_candidate`
- `interchromosomal_rearrangement`

## Cluster-level event types

Breakdance currently emits:

- `complex_rearrangement_cluster`
- `chromothripsis_like_cluster`

These are heuristic labels intended for ranking and review.

---

# Scoring logic

Each interpreted event receives:

- `priority_score`
- `priority_tier`

Scoring currently considers:

- event type
- number of supporting callers
- median support
- confidence level
- cluster complexity
- interchromosomal complexity
- adaptive hotspot hints
- optional focus genes

Priority tiers:

- `Tier1` = highest priority
- `Tier2` = medium priority
- `Tier3` = lower priority / exploratory

---

# Adaptive sequencing support

Breakdance includes a **target-free adaptive-sequencing-aware heuristic**.

If no loci list or target BED is available, Breakdance cannot definitively label
breakpoints as on-target or off-target.

Instead, it computes a local breakpoint-density hint that can help flag:

- hotspot-like breakpoint regions
- targeted-looking breakpoint concentrations
- enriched local rearrangement neighborhoods

This is a **heuristic only**, not a formal targeting label.

---

# Troubleshooting

## `run_harmonize.py` cannot find VCFs
Check:

- BAM-derived sample naming
- whether you used `--output-root`
- whether caller outputs are present under the expected directory
- presence of:
  - `sniffles2/sniffles2.vcf(.gz)`
  - `cutesv/cutesv.vcf(.gz)`
  - `svim/variants.vcf(.gz)`

## BAM list file not found
Check:
- path passed to `--bam-list`
- working directory
- file permissions

## Annotation appears empty
Check:
- BED files are valid BED4
- chromosome naming matches between VCFs and BEDs
  - e.g. `chr8` vs `8`

## Too many separate clusters
Try increasing:
```bash
--cluster-distance-bp
```

## Too many harmonized junctions
Try adjusting:
```bash
--tolerance-bp
```

## `setup_references.py` fails for T2T
This usually means no default source was available.
Provide explicit:
- `--fasta` or `--fasta-url`
- `--gtf` or `--gtf-url`

## Enhancer BED not present
That is expected unless you explicitly provided one.
Enhancers are optional.

---

# Example end-to-end workflow

## Step 1: prepare references
```bash
python setup_references.py \
  --build hg38 \
  --reference-root references
```

## Step 2: make sure caller outputs exist
Confirm for each sample:
- Sniffles2 VCF exists
- cuteSV VCF exists
- SVIM VCF exists

## Step 3: run Breakdance
```bash
python run_harmonize.py \
  --bam-list bams_to_run.txt \
  --output-root /path/to/all_sv_calls \
  --genes-bed references/hg38/annotation/bed/genes.bed \
  --exons-bed references/hg38/annotation/bed/exons.bed \
  --promoters-bed references/hg38/annotation/bed/promoters.bed \
  --summary-dir harmonized_summary
```

## Step 4: inspect outputs
Start with:
- `*.interpreted_events.tsv`
- `*.junction_clusters.tsv`
- `harmonized_cohort_summary.tsv`

## Step 5: prioritize candidates
Review first:
- `Tier1` and `Tier2` events
- multi-caller supported junctions
- large clusters
- interchromosomal clusters
- focus-gene hits

---

# Intended use cases

Breakdance is especially useful for:

- ONT leukemia / lymphoma SV review
- long-read cancer rearrangement triage
- harmonizing multiple ONT SV callers
- identifying candidate fusions
- identifying candidate complex SV loci
- generating structured outputs for downstream graph analysis

---

# Limitations

Breakdance Version 1 does **not currently include**:

- copy-number integration
- CN-balanced graph inference
- segment graph reconstruction
- formal chromothripsis criteria
- formal BFB detection
- formal ecDNA detection
- haplotype-aware breakpoint chaining
- read-level validation
- clinical classification logic

---

# Recommended future extensions

Useful next modules include:

- `copy_number.py`
- `complex_events.py`
- `graph_analysis.py`
- `gene_sets.py`
- centralized config support (`config.yaml`)

---

# Status

Breakdance is currently a **starter research pipeline** and should be treated as
prototype / development code unless you add your own validation and QC controls.

---
