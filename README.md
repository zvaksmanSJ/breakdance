# Breakdance

**Breakdance** is a starter ONT structural variant interpretation pipeline for long-read structural variant callsets.

It harmonizes structural variant VCFs from **Sniffles2**, **cuteSV**, and **SVIM** into richer breakpoint-, event-, fusion-, cluster-, and graph-level outputs.

## Current workflow

For each sample, Breakdance can now:
- parse caller-specific VCF evidence
- harmonize near-equivalent breakpoint junctions across callers
- cluster nearby junctions into candidate event groups
- annotate breakpoint sides and event spans with genes, exons, promoters, and enhancers
- derive fusion candidates and other interpreted SV events
- score and prioritize events and fusions
- export an enriched breakpoint graph
- optionally estimate local breakpoint coverage from BAMs
- optionally create placeholder figure artifacts for downstream review
- write per-sample TSV and JSON outputs plus a cohort summary

## Main modules
- `schema.py` — dataclasses for calls, junctions, clusters, fusions, events, graphs
- `vcf_parser.py` — caller-aware VCF parsing
- `harmonize.py` — cross-caller grouping and consensus junction generation
- `clustering.py` — proximity-based connected-component clustering
- `annotate.py` — breakpoint-side and span-based interval annotations
- `interpret.py` — fusion and structural event interpretation
- `score.py` — event and fusion prioritization
- `graph.py` — enriched breakpoint graph export
- `report.py` — TSV/JSON writers
- `coverage.py` — optional local BAM coverage evidence
- `fusions.py` — fusion consolidation helpers
- `figures.py` — placeholder figure generation
- `lsf.py` — lightweight LSF launcher helper
- `run_harmonize.py` — end-to-end CLI entrypoint
- `setup_references.py` — reference/annotation preparation

## Requirements
Core pipeline uses the Python standard library plus:
- `pysam` for optional BAM-based local coverage estimation

Install:
```bash
pip install -r requirements.txt
```

If `pysam` is unavailable, Breakdance still runs; coverage-related fields remain unpopulated.

## Minimal run
```bash
python run_harmonize.py \
  --bam-list /path/to/bams_to_run.txt \
  --summary-dir harmonized_summary
```

## Run with annotations, focus genes, and figures
```bash
python run_harmonize.py \
  --bam-list /path/to/bams_to_run.txt \
  --output-root /path/to/all_sv_calls \
  --genes-bed references/hg38/annotation/bed/genes.bed \
  --exons-bed references/hg38/annotation/bed/exons.bed \
  --promoters-bed references/hg38/annotation/bed/promoters.bed \
  --enhancers-bed references/hg38/annotation/bed/enhancers.bed \
  --focus-genes MYC,BCL2,TP53 \
  --figure-dir harmonized_summary/figures \
  --summary-dir harmonized_summary
```

## Key outputs per sample
- `<sample>.harmonized_breakpoints.tsv`
- `<sample>.junction_clusters.tsv`
- `<sample>.fusion_events.tsv`
- `<sample>.interpreted_events.tsv`
- `<sample>.breakpoint_graph.json`
- `<sample>.dashboard.json`
- `<sample>.figures_manifest.tsv`

## Notes and limitations
- coverage fields are approximate local depth estimates, not full breakpoint-spanning read validation
- figure generation currently creates lightweight text placeholders for stable downstream references
- enhancer annotations are optional
- this remains a review-oriented research pipeline, not a clinically validated interpretation engine
