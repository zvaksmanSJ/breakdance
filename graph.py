#!/usr/bin/env python3
"""
graph.py

Build a lightweight breakpoint graph from harmonized junctions with enriched
support and event metadata.
"""

from __future__ import annotations

from schema import BreakpointGraph


def build_breakpoint_graph(sample: str, junctions, events=None, fusions=None) -> BreakpointGraph:
    nodes = []
    edges = []
    seen_nodes = set()
    event_by_junction = {}
    fusion_by_junction = {}

    for event in events or []:
        for jid in event.junction_ids:
            event_by_junction.setdefault(jid, []).append(
                {
                    "event_id": event.event_id,
                    "event_type": event.event_type,
                    "event_label": event.event_label,
                    "priority_score": event.priority_score,
                    "priority_tier": event.priority_tier,
                }
            )

    for fusion in fusions or []:
        for jid in fusion.junction_ids:
            fusion_by_junction.setdefault(jid, []).append(
                {
                    "fusion_id": fusion.fusion_id,
                    "fusion_label": fusion.fusion_label,
                    "caller_count": fusion.caller_count,
                    "priority_score": fusion.priority_score,
                    "priority_tier": fusion.priority_tier,
                }
            )

    for j in junctions:
        node1 = f"{j.chrom1}:{j.pos1}"
        node2 = f"{j.chrom2}:{j.pos2}"

        if node1 not in seen_nodes:
            nodes.append(
                {
                    "node_id": node1,
                    "chrom": j.chrom1,
                    "pos": j.pos1,
                    "side": 1,
                    "affected_genes": [g for g in j.affected_genes if g],
                    "affected_exons": list(j.affected_exons),
                }
            )
            seen_nodes.add(node1)

        if node2 not in seen_nodes:
            nodes.append(
                {
                    "node_id": node2,
                    "chrom": j.chrom2,
                    "pos": j.pos2,
                    "side": 2,
                    "affected_genes": [g for g in j.affected_genes if g],
                    "affected_exons": list(j.affected_exons),
                }
            )
            seen_nodes.add(node2)

        edges.append(
            {
                "edge_id": j.junction_id,
                "source": node1,
                "target": node2,
                "svtype": j.svtype,
                "strand1": j.strand1,
                "strand2": j.strand2,
                "caller_count": j.caller_count,
                "callers": list(j.callers),
                "support_by_caller": dict(j.support_by_caller),
                "confidence": j.confidence,
                "is_interchromosomal": j.is_interchromosomal,
                "supporting_reads": j.supporting_reads,
                "covering_reads": j.covering_reads,
                "non_supporting_reads": j.non_supporting_reads,
                "affected_genes": list(j.affected_genes),
                "affected_exons": list(j.affected_exons),
                "event_region1": j.event_region1,
                "event_region2": j.event_region2,
                "events": event_by_junction.get(j.junction_id, []),
                "fusions": fusion_by_junction.get(j.junction_id, []),
            }
        )

    metadata = {
        "node_count": len(nodes),
        "edge_count": len(edges),
        "interchromosomal_edge_count": sum(1 for e in edges if e["is_interchromosomal"]),
        "event_count": len(events or []),
        "fusion_count": len(fusions or []),
    }

    return BreakpointGraph(
        sample=sample,
        nodes=nodes,
        edges=edges,
        metadata=metadata,
    )
