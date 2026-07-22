#!/usr/bin/env python3
"""
score.py

Purpose
-------
Assign priority scores and tiers to interpreted events.

Why this module exists
----------------------
Interpretation can produce many candidate events. Scoring helps rank them
for review by combining:
- event type
- caller support
- read support
- confidence
- cluster complexity
- optional focus genes of interest
"""

from __future__ import annotations


def score_events(
    events,
    junction_by_id,
    cluster_by_id=None,
    focus_genes: set[str] | None = None,
):
    """
    Add priority_score and priority_tier to interpreted events.

    Parameters
    ----------
    events
        List of InterpretedEvent objects.

    junction_by_id
        Mapping from junction_id -> HarmonizedJunction.

    cluster_by_id
        Optional mapping from cluster_id -> SVCluster.

    focus_genes
        Optional uppercase gene set used to boost event priority.
    """
    focus_genes = focus_genes or set()

    for e in events:
        score = 0.0

        # Base weight by event type.
        event_weights = {
            "fusion_candidate": 5.0,
            "exon_disruption": 4.0,
            "promoter_hijack_candidate": 4.5,
            "enhancer_hijack_candidate": 4.0,
            "interchromosomal_rearrangement": 3.0,
            "complex_rearrangement_cluster": 5.5,
            "chromothripsis_like_cluster": 7.0,
        }
        score += event_weights.get(e.event_type, 1.0)

        # Add evidence from supporting junctions.
        for jid in e.junction_ids:
            j = junction_by_id.get(jid)
            if not j:
                continue

            score += min(j.caller_count, 3) * 1.0

            if j.median_support is not None:
                score += min(float(j.median_support), 10.0) * 0.2

            if j.confidence == "high":
                score += 2.0
            elif j.confidence == "medium":
                score += 1.0

            if j.adaptive_hotspot_cluster:
                score += 0.5

        # Add cluster complexity contribution when relevant.
        if e.cluster_id and cluster_by_id and e.cluster_id in cluster_by_id:
            c = cluster_by_id[e.cluster_id]
            score += min(c.junction_count, 10) * 0.3
            score += min(c.interchromosomal_junction_count, 5) * 0.5

        # Boost if any gene is in the focus list.
        if focus_genes:
            if any(g.upper() in focus_genes for g in e.genes):
                score += 3.0

        e.priority_score = round(score, 3)

        if score >= 12:
            e.priority_tier = "Tier1"
        elif score >= 8:
            e.priority_tier = "Tier2"
        else:
            e.priority_tier = "Tier3"

    return events
