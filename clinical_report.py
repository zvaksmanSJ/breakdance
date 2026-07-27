#!/usr/bin/env python3
"""
clinical_report.py

Clinical-style exports for Breakdance, including TSV/JSON summaries and a
self-contained HTML report suitable for analyst review.
"""

from __future__ import annotations

import html
import json
from pathlib import Path


def _safe(value):
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ", ".join(str(x) for x in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True)
    return str(value)


def _event_rows(sample: str, events, fusions):
    rows = []
    fusion_by_id = {f.fusion_id: f for f in fusions}

    for event in events:
        fusion = fusion_by_id.get(event.event_id)
        rows.append(
            {
                "sample": sample,
                "record_id": event.event_id,
                "record_type": "fusion" if event.event_type == "fusion_event" else "event",
                "event_type": event.event_type,
                "event_label": event.event_label or "",
                "cluster_id": event.cluster_id or "",
                "genes": event.genes,
                "affected_genes": event.affected_genes,
                "fusion_label": fusion.fusion_label if fusion else "",
                "chrom_location": event.location_summary,
                "junction_ids": event.junction_ids,
                "supporting_callers": event.supporting_callers,
                "support_by_caller": event.support_by_caller,
                "supporting_reads": event.supporting_reads,
                "covering_reads": event.covering_reads,
                "non_supporting_reads": event.non_supporting_reads,
                "confidence": event.confidence,
                "priority_score": event.priority_score,
                "priority_tier": event.priority_tier,
                "figure_paths": event.figure_paths,
                "details": event.details,
            }
        )

    return rows


def write_clinical_summary_tsv(out_path: Path, sample: str, events, fusions) -> None:
    rows = _event_rows(sample, events, fusions)
    header = [
        "sample",
        "record_id",
        "record_type",
        "event_type",
        "event_label",
        "cluster_id",
        "genes",
        "affected_genes",
        "fusion_label",
        "chrom_location",
        "junction_ids",
        "supporting_callers",
        "support_by_caller",
        "supporting_reads",
        "covering_reads",
        "non_supporting_reads",
        "confidence",
        "priority_score",
        "priority_tier",
        "figure_paths",
        "details",
    ]

    with out_path.open("w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(_safe(row.get(col, "")) for col in header) + "\n")


def write_clinical_summary_json(out_path: Path, sample: str, junctions, clusters, events, fusions, graph) -> None:
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
        "clinical_records": _event_rows(sample, events, fusions),
        "junctions": [
            {
                "junction_id": j.junction_id,
                "svtype": j.svtype,
                "chrom1": j.chrom1,
                "pos1": j.pos1,
                "chrom2": j.chrom2,
                "pos2": j.pos2,
                "callers": j.callers,
                "confidence": j.confidence,
                "affected_genes": j.affected_genes,
                "affected_exons": j.affected_exons,
                "supporting_reads": j.supporting_reads,
                "covering_reads": j.covering_reads,
                "non_supporting_reads": j.non_supporting_reads,
            }
            for j in junctions
        ],
        "clusters": [
            {
                "cluster_id": c.cluster_id,
                "junction_count": c.junction_count,
                "chromosomes": c.chromosomes,
                "affected_genes": c.affected_genes,
                "priority_tier": c.priority_tier,
                "priority_score": c.priority_score,
                "figure_path": c.figure_path,
            }
            for c in clusters
        ],
        "fusions": [
            {
                "fusion_id": f.fusion_id,
                "fusion_label": f.fusion_label,
                "location_summary": f.location_summary,
                "supporting_callers": f.supporting_callers,
                "supporting_reads": f.supporting_reads,
                "covering_reads": f.covering_reads,
                "non_supporting_reads": f.non_supporting_reads,
                "priority_tier": f.priority_tier,
                "priority_score": f.priority_score,
                "figure_path": f.figure_path,
            }
            for f in fusions
        ],
        "graph": graph.metadata,
    }

    with out_path.open("w") as out:
        json.dump(payload, out, indent=2, sort_keys=True)


def _html_table(title: str, rows: list[dict], columns: list[str]) -> str:
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(_safe(row.get(col, '')))}</td>" for col in columns)
        body.append(f"<tr>{cells}</tr>")
    if not body:
        body.append(f"<tr><td colspan=\"{len(columns)}\">No records</td></tr>")
    header = "".join(f"<th>{html.escape(col)}</th>" for col in columns)
    return f"<h2>{html.escape(title)}</h2><table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def write_clinical_summary_html(out_path: Path, sample: str, junctions, clusters, events, fusions, graph) -> None:
    event_rows = _event_rows(sample, events, fusions)
    summary_rows = [
        {"metric": "sample", "value": sample},
        {"metric": "junction_count", "value": len(junctions)},
        {"metric": "cluster_count", "value": len(clusters)},
        {"metric": "event_count", "value": len(events)},
        {"metric": "fusion_count", "value": len(fusions)},
        {"metric": "graph_nodes", "value": len(graph.nodes)},
        {"metric": "graph_edges", "value": len(graph.edges)},
    ]
    junction_rows = [
        {
            "junction_id": j.junction_id,
            "svtype": j.svtype,
            "location": f"{j.chrom1}:{j.pos1} -> {j.chrom2}:{j.pos2}",
            "callers": j.callers,
            "confidence": j.confidence,
            "affected_genes": j.affected_genes,
            "affected_exons": j.affected_exons,
            "supporting_reads": j.supporting_reads,
            "covering_reads": j.covering_reads,
            "non_supporting_reads": j.non_supporting_reads,
        }
        for j in junctions
    ]
    cluster_rows = [
        {
            "cluster_id": c.cluster_id,
            "junction_count": c.junction_count,
            "chromosomes": c.chromosomes,
            "affected_genes": c.affected_genes,
            "priority_tier": c.priority_tier,
            "priority_score": c.priority_score,
            "figure_path": c.figure_path,
        }
        for c in clusters
    ]
    fusion_rows = [
        {
            "fusion_id": f.fusion_id,
            "fusion_label": f.fusion_label,
            "location_summary": f.location_summary,
            "supporting_callers": f.supporting_callers,
            "supporting_reads": f.supporting_reads,
            "covering_reads": f.covering_reads,
            "non_supporting_reads": f.non_supporting_reads,
            "confidence": f.confidence,
            "priority_tier": f.priority_tier,
            "priority_score": f.priority_score,
            "figure_path": f.figure_path,
        }
        for f in fusions
    ]
    graph_rows = [{"metric": key, "value": value} for key, value in sorted(graph.metadata.items())]

    html_text = f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Breakdance clinical summary - {html.escape(sample)}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #222; }}
    h1, h2 {{ color: #0b3d91; }}
    table {{ border-collapse: collapse; width: 100%; margin-bottom: 28px; table-layout: fixed; }}
    th, td {{ border: 1px solid #ccc; padding: 6px 8px; vertical-align: top; word-wrap: break-word; }}
    th {{ background: #f3f6fb; text-align: left; }}
    .note {{ background: #fff8e1; border: 1px solid #e0c97f; padding: 10px; margin-bottom: 24px; }}
  </style>
</head>
<body>
  <h1>Breakdance clinical summary</h1>
  <p><strong>Sample:</strong> {html.escape(sample)}</p>
  <div class="note">
    <strong>Interpretation note:</strong>
    This report is intended for research and structured review support. Coverage values are approximate local estimates unless separately validated.
  </div>
  {_html_table("Sample summary", summary_rows, ["metric", "value"])}
  {_html_table("Clinical event summary", event_rows, ["record_id", "record_type", "event_type", "event_label", "fusion_label", "chrom_location", "genes", "affected_genes", "supporting_callers", "supporting_reads", "covering_reads", "non_supporting_reads", "confidence", "priority_tier", "priority_score", "figure_paths"])}
  {_html_table("Fusion calls", fusion_rows, ["fusion_id", "fusion_label", "location_summary", "supporting_callers", "supporting_reads", "covering_reads", "non_supporting_reads", "confidence", "priority_tier", "priority_score", "figure_path"])}
  {_html_table("Harmonized junctions", junction_rows, ["junction_id", "svtype", "location", "callers", "confidence", "affected_genes", "affected_exons", "supporting_reads", "covering_reads", "non_supporting_reads"])}
  {_html_table("Clusters", cluster_rows, ["cluster_id", "junction_count", "chromosomes", "affected_genes", "priority_tier", "priority_score", "figure_path"])}
  {_html_table("Graph summary", graph_rows, ["metric", "value"])}
</body>
</html>
"""
    with out_path.open("w") as out:
        out.write(html_text)
