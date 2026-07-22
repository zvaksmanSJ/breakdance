#!/usr/bin/env python3
"""
run_harmonize.py

Purpose
-------
Main CLI entrypoint for the Breakdance Version 1 ONT SV interpretation pipeline.

What this script does
---------------------
For each sample identified from a BAM list, this script:

1. Locates caller output VCFs from:
   - Sniffles2
   - cuteSV
   - SVIM

2. Parses raw structural variant calls

3. Adds adaptive-sequencing-aware density hints
   (without requiring a target BED)

4. Harmonizes raw calls into consensus breakpoint junctions

5. Clusters harmonized junctions into candidate event groups

6. Optionally annotates junctions and clusters using BED4 interval files:
   - genes
   - exons
   - promoters
   - enhancers

7. Interprets biological / structural candidate events

8. Scores and prioritizes interpreted events

9. Builds a lightweight breakpoint graph

10. Writes output files per sample

11. Writes a cohort-level summary TSV

Input model
-----------
This follows the same BAM-list style used by your earlier SV caller workflow.

If --output-root is provided:
    <output-root>/<sample>/

Otherwise:
    <bam_parent>/<sample>_sv_callers/

Outputs per sample
------------------
- <sample>.harmonized_breakpoints.tsv
- <sample>.junction_clusters.tsv
- <sample>.interpreted_events.tsv
- <sample>.breakpoint_graph.json
- <sample>.dashboard.json

Cohort summary
--------------
- harmonized_cohort_summary.tsv

Notes
-----
This is a Version 1 interpretation layer:
- useful for prioritization and review
- not yet a full CN-balanced graph reconstruction system
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from annotate import annotate_clusters, annotate_junctions
from clustering import cluster_junctions
from graph import build_breakpoint_graph
from harmonize import annotate_density_hints, harmonize_calls
from interpret import interpret_clusters, interpret_junctions
from report import (
    write_clusters_tsv,
    write_dashboard_json,
    write_events_tsv,
    write_graph_json,
    write_harmonized_tsv,
)
from score import score_events
from vcf_parser import discover_caller_vcfs, read_vcf_calls


def log(msg: str) -> None:
    """
    Print a progress message immediately.
    """
    print(msg, flush=True)


def fail(msg: str, code: int = 1) -> None:
    """
    Print an error message and exit.
    """
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


def sample_name_from_bam(bam_path: Path) -> str:
    """
    Derive a sample name from a BAM filename.

    Example
    -------
    sample1.bam -> sample1
    """
    name = bam_path.name
    if name.endswith(".bam"):
        return name[:-4]
    return bam_path.stem


def load_bam_list(bam_list_path: Path) -> list[Path]:
    """
    Load BAM file paths from a text file.

    Rules
    -----
    - one BAM path per line
    - blank lines ignored
    - comment lines starting with '#' ignored
    - listed BAMs must exist
    """
    if not bam_list_path.exists():
        fail(f"BAM list file not found: {bam_list_path}")

    bam_paths: list[Path] = []

    with bam_list_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            bam = Path(line)

            if not bam.exists():
                fail(f"BAM from list does not exist: {bam}")
            if not bam.is_file():
                fail(f"BAM from list is not a file: {bam}")

            bam_paths.append(bam)

    if not bam_paths:
        fail(f"No BAMs found in list file: {bam_list_path}")

    return bam_paths


def resolve_sample_output_dir(
    bam_path: Path,
    sample: str,
    output_root: str | None,
) -> Path:
    """
    Resolve the sample output directory where caller VCFs are expected.

    If output_root is provided:
        <output_root>/<sample>

    Else:
        <bam_parent>/<sample>_sv_callers
    """
    if output_root:
        return Path(output_root) / sample

    return bam_path.parent / f"{sample}_sv_callers"


def write_summary_tsv(out_path: Path, rows: list[dict[str, str]]) -> None:
    """
    Write a cohort-level summary TSV.
    """
    header = [
        "sample",
        "bam",
        "sample_output_dir",
        "vcf_callers_found",
        "raw_call_count",
        "harmonized_junction_count",
        "cluster_count",
        "event_count",
        "largest_cluster_junction_count",
        "harmonized_tsv",
        "clusters_tsv",
        "events_tsv",
        "graph_json",
        "dashboard_json",
    ]

    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row.get(col, "") for col in header) + "\n")


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Breakdance Version 1 ONT SV interpretation pipeline."
    )

    parser.add_argument(
        "--bam-list",
        required=True,
        help="Text file with one BAM path per line",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="Optional common output root used by the caller workflow",
    )
    parser.add_argument(
        "--tolerance-bp",
        type=int,
        default=100,
        help="Breakpoint tolerance for harmonization",
    )
    parser.add_argument(
        "--cluster-distance-bp",
        type=int,
        default=10000,
        help="Breakpoint distance threshold for clustering harmonized junctions",
    )
    parser.add_argument(
        "--summary-dir",
        default=".",
        help="Directory for cohort summary output",
    )

    parser.add_argument(
        "--genes-bed",
        default=None,
        help="Optional BED4 genes annotation file",
    )
    parser.add_argument(
        "--exons-bed",
        default=None,
        help="Optional BED4 exons annotation file",
    )
    parser.add_argument(
        "--promoters-bed",
        default=None,
        help="Optional BED4 promoters annotation file",
    )
    parser.add_argument(
        "--enhancers-bed",
        default=None,
        help="Optional BED4 enhancers annotation file",
    )
    parser.add_argument(
        "--focus-genes",
        default="",
        help="Comma-separated gene list used to boost event ranking",
    )

    return parser.parse_args()


def main() -> None:
    """
    Run the full Breakdance Version 1 pipeline.
    """
    args = parse_args()

    bam_paths = load_bam_list(Path(args.bam_list))
    summary_rows: list[dict[str, str]] = []

    focus_genes = set()
    if args.focus_genes.strip():
        focus_genes = {
            x.strip().upper()
            for x in args.focus_genes.split(",")
            if x.strip()
        }

    log("============================================================")
    log("Breakdance Version 1 ONT SV interpretation pipeline")
    log("============================================================")
    log(f"BAM list              : {args.bam_list}")
    log(f"Output root           : {args.output_root if args.output_root else '[same directory as each BAM]'}")
    log(f"Harmonize tolerance   : {args.tolerance_bp}")
    log(f"Cluster distance (bp) : {args.cluster_distance_bp}")
    log(f"Summary dir           : {args.summary_dir}")
    log(f"Genes BED             : {args.genes_bed}")
    log(f"Exons BED             : {args.exons_bed}")
    log(f"Promoters BED         : {args.promoters_bed}")
    log(f"Enhancers BED         : {args.enhancers_bed}")
    log(f"Focus genes           : {','.join(sorted(focus_genes)) if focus_genes else '[none]'}")
    log("============================================================")

    for bam in bam_paths:
        sample = sample_name_from_bam(bam)
        sample_output_dir = resolve_sample_output_dir(
            bam_path=bam,
            sample=sample,
            output_root=args.output_root,
        )

        found_vcfs = discover_caller_vcfs(sample_output_dir)

        log("------------------------------------------------------------")
        log(f"Sample            : {sample}")
        log(f"Input BAM         : {bam}")
        log(f"Sample output dir : {sample_output_dir}")
        log(f"VCFs found        : {', '.join(sorted(found_vcfs)) if found_vcfs else 'none'}")

        raw_calls = []

        for caller, vcf_path in found_vcfs.items():
            calls = read_vcf_calls(vcf_path, caller, sample)
            log(f"  {caller}: parsed {len(calls)} calls from {vcf_path}")
            raw_calls.extend(calls)

        if not raw_calls:
            summary_rows.append(
                {
                    "sample": sample,
                    "bam": str(bam),
                    "sample_output_dir": str(sample_output_dir),
                    "vcf_callers_found": "",
                    "raw_call_count": "0",
                    "harmonized_junction_count": "0",
                    "cluster_count": "0",
                    "event_count": "0",
                    "largest_cluster_junction_count": "0",
                    "harmonized_tsv": "",
                    "clusters_tsv": "",
                    "events_tsv": "",
                    "graph_json": "",
                    "dashboard_json": "",
                }
            )
            continue

        # Step 1: adaptive-density hints
        raw_calls = annotate_density_hints(raw_calls)

        # Step 2: harmonization
        junctions = harmonize_calls(
            raw_calls,
            tolerance_bp=args.tolerance_bp,
        )

        # Step 3: clustering
        clusters = cluster_junctions(
            junctions,
            distance_bp=args.cluster_distance_bp,
        )

        # Step 4: interval annotation
        junctions = annotate_junctions(
            junctions,
            genes_bed=args.genes_bed,
            exons_bed=args.exons_bed,
            promoters_bed=args.promoters_bed,
            enhancers_bed=args.enhancers_bed,
        )

        junction_by_id = {j.junction_id: j for j in junctions}

        clusters = annotate_clusters(
            clusters,
            junction_by_id=junction_by_id,
        )

        cluster_by_id = {c.cluster_id: c for c in clusters}

        # Step 5: interpretation
        events = (
            interpret_junctions(sample, junctions)
            + interpret_clusters(sample, clusters)
        )

        # Step 6: scoring
        events = score_events(
            events,
            junction_by_id=junction_by_id,
            cluster_by_id=cluster_by_id,
            focus_genes=focus_genes,
        )

        # Sort higher-priority events first
        events.sort(
            key=lambda e: (-e.priority_score, e.event_type, e.event_id)
        )

        # Step 7: graph export
        graph = build_breakpoint_graph(sample, junctions)

        # Output paths
        harmonized_tsv = sample_output_dir / f"{sample}.harmonized_breakpoints.tsv"
        clusters_tsv = sample_output_dir / f"{sample}.junction_clusters.tsv"
        events_tsv = sample_output_dir / f"{sample}.interpreted_events.tsv"
        graph_json = sample_output_dir / f"{sample}.breakpoint_graph.json"
        dashboard_json = sample_output_dir / f"{sample}.dashboard.json"

        # Step 8: write outputs
        write_harmonized_tsv(harmonized_tsv, junctions)
        write_clusters_tsv(clusters_tsv, clusters)
        write_events_tsv(events_tsv, events)
        write_graph_json(graph_json, graph)
        write_dashboard_json(dashboard_json, sample, junctions, clusters, events, graph)

        log(f"Wrote harmonized TSV : {harmonized_tsv}")
        log(f"Wrote clusters TSV   : {clusters_tsv}")
        log(f"Wrote events TSV     : {events_tsv}")
        log(f"Wrote graph JSON     : {graph_json}")
        log(f"Wrote dashboard JSON : {dashboard_json}")
        log(f"Raw calls            : {len(raw_calls)}")
        log(f"Harmonized junctions : {len(junctions)}")
        log(f"Clusters             : {len(clusters)}")
        log(f"Events               : {len(events)}")

        largest_cluster = max((c.junction_count for c in clusters), default=0)

        summary_rows.append(
            {
                "sample": sample,
                "bam": str(bam),
                "sample_output_dir": str(sample_output_dir),
                "vcf_callers_found": ",".join(sorted(found_vcfs)),
                "raw_call_count": str(len(raw_calls)),
                "harmonized_junction_count": str(len(junctions)),
                "cluster_count": str(len(clusters)),
                "event_count": str(len(events)),
                "largest_cluster_junction_count": str(largest_cluster),
                "harmonized_tsv": str(harmonized_tsv),
                "clusters_tsv": str(clusters_tsv),
                "events_tsv": str(events_tsv),
                "graph_json": str(graph_json),
                "dashboard_json": str(dashboard_json),
            }
        )

    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = summary_dir / "harmonized_cohort_summary.tsv"
    write_summary_tsv(summary_tsv, summary_rows)

    log("============================================================")
    log(f"Wrote cohort summary: {summary_tsv}")
    log("Done.")
    log("============================================================")


if __name__ == "__main__":
    main()
