#!/usr/bin/env python3
"""
report.py

Purpose
-------
Write standardized TSV and JSON outputs from the Breakdance pipeline.

Why this module exists
----------------------
Keeping file-writing logic in one place makes output formats:
- easier to maintain
- easier to version
- easier to reuse in dashboard and downstream analysis
"""

from __future__ import annotations

import json
from pathlib import Path


def _jsonish(value):
    """
    Convert dict/list values into compact JSON strings for TSV output.
    """
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return value


def write_harmonized_tsv(out_path: Path, junctions) -> None:
    """
    Write harmonized junctions to TSV.
    """
    header = [
        "junction_id",
        "svtype",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "strand1",
        "strand2",
        "callers",
        "caller_count",
        "median_support",
        "confidence",
        "chrom_pair",
        "is_interchromosomal",
        "adaptive_hotspot_cluster",
        "adaptive_density_score",
        "raw_call_ids",
        "support_by_caller",
        "qual_by_caller",
        "filters_by_caller",
        "annotations",
        "source_files",
    ]

    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")

        for j in junctions:
            row = [
                j.junction_id,
                j.svtype,
                j.chrom1,
                str(j.pos1),
                j.chrom2,
                str(j.pos2),
                j.strand1 or "",
                j.strand2 or "",
                ",".join(j.callers),
                str(j.caller_count),
                "" if j.median_support is None else str(j.median_support),
                j.confidence,
                j.chrom_pair or "",
                str(j.is_interchromosomal),
                str(j.adaptive_hotspot_cluster),
                "" if j.adaptive_density_score is None else str(j.adaptive_density_score),
                ",".join(j.raw_call_ids),
                _jsonish(j.support_by_caller),
                _jsonish(j.qual_by_caller),
                _jsonish(j.filters_by_caller),
                _jsonish(j.annotations),
                ";".join(j.source_files),
            ]
            out.write("\t".join(map(str, row)) + "\n")


def write_clusters_tsv(out_path: Path, clusters) -> None:
    """
    Write cluster summaries to TSV.
    """
    header = [
        "cluster_id",
        "junction_ids",
        "junction_count",
        "chromosomes",
        "interchromosomal_junction_count",
        "caller_set",
        "max_adaptive_density_score",
        "any_adaptive_hotspot",
        "min_position_by_chrom",
        "max_position_by_chrom",
        "notes",
        "annotations",
    ]

    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")

        for c in clusters:
            row = [
                c.cluster_id,
                ",".join(c.junction_ids),
                str(c.junction_count),
                ",".join(c.chromosomes),
                str(c.interchromosomal_junction_count),
                ",".join(c.caller_set),
                "" if c.max_adaptive_density_score is None else str(c.max_adaptive_density_score),
                str(c.any_adaptive_hotspot),
                _jsonish(c.min_position_by_chrom),
                _jsonish(c.max_position_by_chrom),
                _jsonish(c.notes),
                _jsonish(c.annotations),
            ]
            out.write("\t".join(map(str, row)) + "\n")


def write_events_tsv(out_path: Path, events) -> None:
    """
    Write interpreted events to TSV.
    """
    header = [
        "event_id",
        "event_type",
        "sample",
        "cluster_id",
        "junction_ids",
        "genes",
        "confidence",
        "priority_score",
        "priority_tier",
        "details",
    ]

    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")

        for e in events:
            row = [
                e.event_id,
                e.event_type,
                e.sample,
                e.cluster_id or "",
                ",".join(e.junction_ids),
                ",".join(e.genes),
                e.confidence,
                str(e.priority_score),
                e.priority_tier,
                _jsonish(e.details),
            ]
            out.write("\t".join(map(str, row)) + "\n")


def write_graph_json(out_path: Path, graph) -> None:
    """
    Write breakpoint graph to JSON.
    """
    payload = {
        "sample": graph.sample,
        "nodes": graph.nodes,
        "edges": graph.edges,
        "metadata": graph.metadata,
    }

    with out_path.open("w") as out:
        json.dump(payload, out, indent=2, sort_keys=True)


def write_dashboard_json(out_path: Path, sample: str, junctions, clusters, events, graph) -> None:
    """
    Write one JSON bundle for dashboard consumption.
    """
    payload = {
        "sample": sample,
        "summary": {
            "junction_count": len(junctions),
            "cluster_count": len(clusters),
            "event_count": len(events),
            "graph_nodes": len(graph.nodes),
            "graph_edges": len(graph.edges),
        },
        "junctions": [
            {
                "junction_id": j.junction_id,
                "svtype": j.svtype,
                "chrom1": j.chrom1,
                "pos1": j.pos1,
                "chrom2": j.chrom2,
                "pos2": j.pos2,
                "callers": j.callers,
                "caller_count": j.caller_count,
                "confidence": j.confidence,
                "annotations": j.annotations,
            }
            for j in junctions
        ],
        "clusters": [
            {
                "cluster_id": c.cluster_id,
                "junction_ids": c.junction_ids,
                "junction_count": c.junction_count,
                "chromosomes": c.chromosomes,
                "notes": c.notes,
                "annotations": c.annotations,
            }
            for c in clusters
        ],
        "events": [
            {
                "event_id": e.event_id,
                "event_type": e.event_type,
                "cluster_id": e.cluster_id,
                "junction_ids": e.junction_ids,
                "genes": e.genes,
                "confidence": e.confidence,
                "priority_score": e.priority_score,
                "priority_tier": e.priority_tier,
                "details": e.details,
            }
            for e in events
        ],
        "graph": {
            "nodes": graph.nodes,
            "edges": graph.edges,
            "metadata": graph.metadata,
        },
    }

    with out_path.open("w") as out:
        json.dump(payload, out, indent=2, sort_keys=True)
