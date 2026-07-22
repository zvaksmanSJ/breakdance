#!/usr/bin/env python3
"""
report.py

Write standardized TSV and JSON outputs from the Breakdance pipeline.
"""

from __future__ import annotations

import json
from pathlib import Path


def _jsonish(value):
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return value


def write_harmonized_tsv(out_path: Path, junctions) -> None:
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
        "consensus_event_type",
        "event_region1",
        "event_region2",
        "affected_genes",
        "affected_exons",
        "affected_promoters",
        "affected_enhancers",
        "supporting_reads",
        "covering_reads",
        "non_supporting_reads",
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
                j.consensus_event_type or "",
                j.event_region1 or "",
                j.event_region2 or "",
                ",".join(j.affected_genes),
                ",".join(j.affected_exons),
                ",".join(j.affected_promoters),
                ",".join(j.affected_enhancers),
                "" if j.supporting_reads is None else str(j.supporting_reads),
                "" if j.covering_reads is None else str(j.covering_reads),
                "" if j.non_supporting_reads is None else str(j.non_supporting_reads),
                ",".join(j.raw_call_ids),
                _jsonish(j.support_by_caller),
                _jsonish(j.qual_by_caller),
                _jsonish(j.filters_by_caller),
                _jsonish(j.annotations),
                ";".join(j.source_files),
            ]
            out.write("\t".join(map(str, row)) + "\n")


def write_clusters_tsv(out_path: Path, clusters) -> None:
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
        "event_span_by_chrom",
        "affected_genes",
        "affected_exons",
        "affected_promoters",
        "affected_enhancers",
        "figure_path",
        "priority_score",
        "priority_tier",
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
                _jsonish(c.event_span_by_chrom),
                ",".join(c.affected_genes),
                ",".join(c.affected_exons),
                ",".join(c.affected_promoters),
                ",".join(c.affected_enhancers),
                c.figure_path or "",
                str(c.priority_score),
                c.priority_tier,
                _jsonish(c.notes),
                _jsonish(c.annotations),
            ]
            out.write("\t".join(map(str, row)) + "\n")


def write_fusions_tsv(out_path: Path, fusions) -> None:
    header = [
        "fusion_id",
        "sample",
        "fusion_label",
        "gene1",
        "gene2",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "orientation",
        "svtype",
        "junction_ids",
        "supporting_callers",
        "support_by_caller",
        "caller_count",
        "supporting_reads",
        "covering_reads",
        "non_supporting_reads",
        "affected_exons",
        "affected_genes",
        "location_summary",
        "figure_path",
        "confidence",
        "priority_score",
        "priority_tier",
        "details",
    ]

    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")
        for f in fusions:
            row = [
                f.fusion_id,
                f.sample,
                f.fusion_label,
                f.gene1,
                f.gene2,
                f.chrom1,
                str(f.pos1),
                f.chrom2,
                str(f.pos2),
                f.orientation or "",
                f.svtype,
                ",".join(f.junction_ids),
                ",".join(f.supporting_callers),
                _jsonish(f.support_by_caller),
                str(f.caller_count),
                "" if f.supporting_reads is None else str(f.supporting_reads),
                "" if f.covering_reads is None else str(f.covering_reads),
                "" if f.non_supporting_reads is None else str(f.non_supporting_reads),
                ",".join(f.affected_exons),
                ",".join(f.affected_genes),
                f.location_summary,
                f.figure_path or "",
                f.confidence,
                str(f.priority_score),
                f.priority_tier,
                _jsonish(f.details),
            ]
            out.write("\t".join(map(str, row)) + "\n")


def write_events_tsv(out_path: Path, events) -> None:
    header = [
        "event_id",
        "event_type",
        "sample",
        "event_label",
        "cluster_id",
        "junction_ids",
        "location_summary",
        "genes",
        "affected_genes",
        "supporting_callers",
        "support_by_caller",
        "supporting_reads",
        "covering_reads",
        "non_supporting_reads",
        "figure_paths",
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
                e.event_label or "",
                e.cluster_id or "",
                ",".join(e.junction_ids),
                e.location_summary,
                ",".join(e.genes),
                ",".join(e.affected_genes),
                ",".join(e.supporting_callers),
                _jsonish(e.support_by_caller),
                "" if e.supporting_reads is None else str(e.supporting_reads),
                "" if e.covering_reads is None else str(e.covering_reads),
                "" if e.non_supporting_reads is None else str(e.non_supporting_reads),
                ";".join(e.figure_paths),
                e.confidence,
                str(e.priority_score),
                e.priority_tier,
                _jsonish(e.details),
            ]
            out.write("\t".join(map(str, row)) + "\n")


def write_graph_json(out_path: Path, graph) -> None:
    payload = {
        "sample": graph.sample,
        "nodes": graph.nodes,
        "edges": graph.edges,
        "metadata": graph.metadata,
    }
    with out_path.open("w") as out:
        json.dump(payload, out, indent=2, sort_keys=True)


def write_figures_manifest_tsv(out_path: Path, sample: str, events, fusions) -> None:
    header = ["sample", "record_type", "record_id", "label", "figure_path"]
    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")
        for e in events:
            for path in e.figure_paths:
                out.write("\t".join([sample, "event", e.event_id, e.event_label or e.event_type, path]) + "\n")
        for f in fusions:
            if f.figure_path:
                out.write("\t".join([sample, "fusion", f.fusion_id, f.fusion_label, f.figure_path]) + "\n")


def write_dashboard_json(out_path: Path, sample: str, junctions, clusters, events, graph, fusions=None) -> None:
    fusions = fusions or []
    payload = {
        "sample": sample,
        "summary": {
            "junction_count": len(junctions),
            "cluster_count": len(clusters),
            "event_count": len(events),
            "fusion_count": len(fusions),
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
                "affected_genes": j.affected_genes,
                "affected_exons": j.affected_exons,
                "supporting_reads": j.supporting_reads,
                "covering_reads": j.covering_reads,
                "non_supporting_reads": j.non_supporting_reads,
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
                "event_span_by_chrom": c.event_span_by_chrom,
                "affected_genes": c.affected_genes,
                "notes": c.notes,
                "annotations": c.annotations,
            }
            for c in clusters
        ],
        "fusions": [
            {
                "fusion_id": f.fusion_id,
                "fusion_label": f.fusion_label,
                "gene1": f.gene1,
                "gene2": f.gene2,
                "location_summary": f.location_summary,
                "supporting_callers": f.supporting_callers,
                "support_by_caller": f.support_by_caller,
                "supporting_reads": f.supporting_reads,
                "affected_exons": f.affected_exons,
                "priority_score": f.priority_score,
                "priority_tier": f.priority_tier,
                "figure_path": f.figure_path,
            }
            for f in fusions
        ],
        "events": [
            {
                "event_id": e.event_id,
                "event_type": e.event_type,
                "event_label": e.event_label,
                "cluster_id": e.cluster_id,
                "junction_ids": e.junction_ids,
                "location_summary": e.location_summary,
                "genes": e.genes,
                "affected_genes": e.affected_genes,
                "supporting_callers": e.supporting_callers,
                "support_by_caller": e.support_by_caller,
                "supporting_reads": e.supporting_reads,
                "priority_score": e.priority_score,
                "priority_tier": e.priority_tier,
                "figure_paths": e.figure_paths,
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
