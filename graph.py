#!/usr/bin/env python3
"""
graph.py

Purpose
-------
Build a lightweight breakpoint graph from harmonized junctions.

Version 1 graph model
---------------------
- node: one breakpoint endpoint
- edge: one harmonized junction connecting two endpoints

Why this module exists
----------------------
This provides a structured representation that can be:
- exported as JSON
- inspected later
- upgraded into a richer segment/junction graph in future versions
"""

from __future__ import annotations

from schema import BreakpointGraph


def build_breakpoint_graph(sample: str, junctions) -> BreakpointGraph:
    """
    Build a breakpoint graph object from harmonized junctions.

    Returns
    -------
    BreakpointGraph
    """
    nodes = []
    edges = []
    seen_nodes = set()

    for j in junctions:
        node1 = f"{j.chrom1}:{j.pos1}"
        node2 = f"{j.chrom2}:{j.pos2}"

        if node1 not in seen_nodes:
            nodes.append(
                {
                    "node_id": node1,
                    "chrom": j.chrom1,
                    "pos": j.pos1,
                }
            )
            seen_nodes.add(node1)

        if node2 not in seen_nodes:
            nodes.append(
                {
                    "node_id": node2,
                    "chrom": j.chrom2,
                    "pos": j.pos2,
                }
            )
            seen_nodes.add(node2)

        edges.append(
            {
                "edge_id": j.junction_id,
                "source": node1,
                "target": node2,
                "svtype": j.svtype,
                "caller_count": j.caller_count,
                "confidence": j.confidence,
                "is_interchromosomal": j.is_interchromosomal,
            }
        )

    return BreakpointGraph(
        sample=sample,
        nodes=nodes,
        edges=edges,
        metadata={
            "node_count": len(nodes),
            "edge_count": len(edges),
        },
    )
