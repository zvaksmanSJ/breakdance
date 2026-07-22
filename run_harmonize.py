#!/usr/bin/env python3
"""
run_harmonize.py

Purpose
-------
Main CLI entrypoint for the Breakdance Version 1 ONT SV interpretation pipeline.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from annotate import annotate_clusters, annotate_junctions
from clustering import cluster_junctions
from coverage import attach_read_evidence, propagate_read_evidence_to_events, propagate_read_evidence_to_fusions
from figures import attach_cluster_figures, attach_event_figures, attach_fusion_figures
from fusions import consolidate_fusions
from graph import build_breakpoint_graph
from harmonize import annotate_density_hints, harmonize_calls
from interpret import interpret_clusters, interpret_junctions
from report import (
    write_clusters_tsv,
    write_dashboard_json,
    write_events_tsv,
    write_figures_manifest_tsv,
    write_fusions_tsv,
    write_graph_json,
    write_harmonized_tsv,
)
from score import score_events, score_fusions
from vcf_parser import discover_caller_vcfs, read_vcf_calls


def log(msg: str) -> None:
    print(msg, flush=True)


def fail(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


def sample_name_from_bam(bam_path: Path) -> str:
    name = bam_path.name
    if name.endswith(".bam"):
        return name[:-4]
    return bam_path.stem


def load_bam_list(bam_list_path: Path) -> list[Path]:
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


def resolve_sample_output_dir(bam_path: Path, sample: str, output_root: str | None) -> Path:
    if output_root:
        return Path(output_root) / sample
    return bam_path.parent / f"{sample}_sv_callers"


def write_summary_tsv(out_path: Path, rows: list[dict[str, str]]) -> None:
    header = [
        "sample",
        "bam",
        "sample_output_dir",
        "vcf_callers_found",
        "raw_call_count",
        "harmonized_junction_count",
        "cluster_count",
        "fusion_count",
        "event_count",
        "largest_cluster_junction_count",
        "harmonized_tsv",
        "clusters_tsv",
        "fusions_tsv",
        "events_tsv",
        "graph_json",
        "dashboard_json",
        "figures_manifest_tsv",
    ]

    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row.get(col, "") for col in header) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Breakdance Version 1 ONT SV interpretation pipeline.")
    parser.add_argument("--bam-list", required=True, help="Text file with one BAM path per line")
    parser.add_argument("--output-root", default=None, help="Optional common output root used by the caller workflow")
    parser.add_argument("--tolerance-bp", type=int, default=100, help="Breakpoint tolerance for harmonization")
    parser.add_argument("--cluster-distance-bp", type=int, default=10000, help="Breakpoint distance threshold for clustering harmonized junctions")
    parser.add_argument("--summary-dir", default=".", help="Directory for cohort summary output")
    parser.add_argument("--genes-bed", default=None, help="Optional BED4 genes annotation file")
    parser.add_argument("--exons-bed", default=None, help="Optional BED4 exons annotation file")
    parser.add_argument("--promoters-bed", default=None, help="Optional BED4 promoters annotation file")
    parser.add_argument("--enhancers-bed", default=None, help="Optional BED4 enhancers annotation file")
    parser.add_argument("--focus-genes", default="", help="Comma-separated gene list used to boost event ranking")
    parser.add_argument("--figure-dir", default=None, help="Optional directory for placeholder event/fusion/cluster figures")
    parser.add_argument("--coverage-window-bp", type=int, default=250, help="Window used for optional local coverage estimation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bam_paths = load_bam_list(Path(args.bam_list))
    summary_rows: list[dict[str, str]] = []

    focus_genes = set()
    if args.focus_genes.strip():
        focus_genes = {x.strip().upper() for x in args.focus_genes.split(",") if x.strip()}

    for bam in bam_paths:
        sample = sample_name_from_bam(bam)
        sample_output_dir = resolve_sample_output_dir(bam_path=bam, sample=sample, output_root=args.output_root)
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        found_vcfs = discover_caller_vcfs(sample_output_dir)

        raw_calls = []
        for caller, vcf_path in found_vcfs.items():
            calls = read_vcf_calls(vcf_path, caller, sample)
            log(f"{sample}: {caller} parsed {len(calls)} calls from {vcf_path}")
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
                    "fusion_count": "0",
                    "event_count": "0",
                    "largest_cluster_junction_count": "0",
                    "harmonized_tsv": "",
                    "clusters_tsv": "",
                    "fusions_tsv": "",
                    "events_tsv": "",
                    "graph_json": "",
                    "dashboard_json": "",
                    "figures_manifest_tsv": "",
                }
            )
            continue

        raw_calls = annotate_density_hints(raw_calls)
        junctions = harmonize_calls(raw_calls, tolerance_bp=args.tolerance_bp)
        junctions = attach_read_evidence(junctions, bam, window_bp=args.coverage_window_bp)

        clusters = cluster_junctions(junctions, distance_bp=args.cluster_distance_bp)
        junctions = annotate_junctions(
            junctions,
            genes_bed=args.genes_bed,
            exons_bed=args.exons_bed,
            promoters_bed=args.promoters_bed,
            enhancers_bed=args.enhancers_bed,
        )
        junction_by_id = {j.junction_id: j for j in junctions}
        clusters = annotate_clusters(clusters, junction_by_id=junction_by_id)
        cluster_by_id = {c.cluster_id: c for c in clusters}

        fusion_events, junction_events = interpret_junctions(sample, junctions)
        cluster_events = interpret_clusters(sample, clusters)
        events = junction_events + cluster_events

        events = propagate_read_evidence_to_events(events, junction_by_id)
        fusion_events = propagate_read_evidence_to_fusions(fusion_events, junction_by_id)

        events = score_events(events, junction_by_id=junction_by_id, cluster_by_id=cluster_by_id, focus_genes=focus_genes)
        fusion_events = score_fusions(fusion_events, focus_genes=focus_genes)
        fusion_events = consolidate_fusions(fusion_events)

        events.sort(key=lambda e: (-e.priority_score, e.event_type, e.event_id))
        fusion_events.sort(key=lambda f: (-f.priority_score, f.fusion_label, f.fusion_id))

        clusters = attach_cluster_figures(sample, clusters, args.figure_dir)
        events = attach_event_figures(sample, events, args.figure_dir)
        fusion_events = attach_fusion_figures(sample, fusion_events, args.figure_dir)

        graph = build_breakpoint_graph(sample, junctions, events=events, fusions=fusion_events)

        harmonized_tsv = sample_output_dir / f"{sample}.harmonized_breakpoints.tsv"
        clusters_tsv = sample_output_dir / f"{sample}.junction_clusters.tsv"
        fusions_tsv = sample_output_dir / f"{sample}.fusion_events.tsv"
        events_tsv = sample_output_dir / f"{sample}.interpreted_events.tsv"
        graph_json = sample_output_dir / f"{sample}.breakpoint_graph.json"
        dashboard_json = sample_output_dir / f"{sample}.dashboard.json"
        figures_manifest_tsv = sample_output_dir / f"{sample}.figures_manifest.tsv"

        write_harmonized_tsv(harmonized_tsv, junctions)
        write_clusters_tsv(clusters_tsv, clusters)
        write_fusions_tsv(fusions_tsv, fusion_events)
        write_events_tsv(events_tsv, events)
        write_graph_json(graph_json, graph)
        write_figures_manifest_tsv(figures_manifest_tsv, sample, events, fusion_events)
        write_dashboard_json(dashboard_json, sample, junctions, clusters, events, graph, fusions=fusion_events)

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
                "fusion_count": str(len(fusion_events)),
                "event_count": str(len(events)),
                "largest_cluster_junction_count": str(largest_cluster),
                "harmonized_tsv": str(harmonized_tsv),
                "clusters_tsv": str(clusters_tsv),
                "fusions_tsv": str(fusions_tsv),
                "events_tsv": str(events_tsv),
                "graph_json": str(graph_json),
                "dashboard_json": str(dashboard_json),
                "figures_manifest_tsv": str(figures_manifest_tsv),
            }
        )

    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = summary_dir / "harmonized_cohort_summary.tsv"
    write_summary_tsv(summary_tsv, summary_rows)
    log(f"Wrote cohort summary: {summary_tsv}")


if __name__ == "__main__":
    main()
