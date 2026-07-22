#!/usr/bin/env python3
"""
score.py

Assign priority scores and tiers to interpreted and fusion events.
"""

from __future__ import annotations


def _junction_evidence_score(j) -> float:
    score = 0.0
    score += min(j.caller_count, 3) * 1.5
    if j.median_support is not None:
        score += min(float(j.median_support), 20.0) * 0.2
    if j.supporting_reads is not None:
        score += min(float(j.supporting_reads), 20.0) * 0.25
    if j.covering_reads is not None and j.supporting_reads is not None and j.covering_reads > 0:
        ratio = j.supporting_reads / max(j.covering_reads, 1)
        score += min(ratio, 1.0) * 2.0
    if j.confidence == "high":
        score += 2.5
    elif j.confidence == "medium":
        score += 1.25
    if j.adaptive_hotspot_cluster:
        score += 0.5
    if j.affected_exons:
        score += 4.0
    if j.is_interchromosomal:
        score += 1.0
    return score


def score_events(
    events,
    junction_by_id,
    cluster_by_id=None,
    focus_genes: set[str] | None = None,
):
    focus_genes = focus_genes or set()

    event_weights = {
        "fusion_event": 8.0,
        "exon_interfering_sv": 7.0,
        "promoter_hijack_candidate": 5.5,
        "enhancer_hijack_candidate": 5.0,
        "interchromosomal_rearrangement": 4.5,
        "large_deletion": 5.0,
        "complex_rearrangement_cluster": 6.5,
        "chromothripsis_like_cluster": 8.5,
    }

    for e in events:
        score = event_weights.get(e.event_type, 2.0)

        for jid in e.junction_ids:
            j = junction_by_id.get(jid)
            if not j:
                continue
            score += _junction_evidence_score(j)

        if e.cluster_id and cluster_by_id and e.cluster_id in cluster_by_id:
            c = cluster_by_id[e.cluster_id]
            score += min(c.junction_count, 12) * 0.35
            score += min(c.interchromosomal_junction_count, 6) * 0.75
            if len(c.affected_genes) >= 2:
                score += 1.0
            if len(c.affected_exons) >= 1:
                score += 2.0

        if e.event_type == "fusion_event":
            score += min(len(e.supporting_callers), 3) * 1.5
            if e.supporting_reads is not None:
                score += min(float(e.supporting_reads), 25.0) * 0.2
            if e.details.get("affected_exons"):
                score += 3.0
            if "<->" in e.location_summary:
                score += 1.5

        if e.event_type == "large_deletion":
            score += min(len(e.affected_genes), 10) * 0.35
            if len(e.affected_genes) >= 3:
                score += 1.0

        if e.event_type == "exon_interfering_sv":
            score += 2.5

        if focus_genes:
            all_genes = {g.upper() for g in (e.genes + e.affected_genes)}
            matching = [g for g in all_genes if g in focus_genes]
            score += len(matching) * 2.5

        e.priority_score = round(score, 3)
        if score >= 18:
            e.priority_tier = "Tier1"
        elif score >= 11:
            e.priority_tier = "Tier2"
        else:
            e.priority_tier = "Tier3"

    return events


def score_fusions(fusions, focus_genes: set[str] | None = None):
    focus_genes = focus_genes or set()
    for fusion in fusions:
        score = 8.0
        score += min(fusion.caller_count, 3) * 1.75
        if fusion.supporting_reads is not None:
            score += min(float(fusion.supporting_reads), 25.0) * 0.25
        if fusion.covering_reads is not None and fusion.supporting_reads is not None and fusion.covering_reads > 0:
            ratio = fusion.supporting_reads / max(fusion.covering_reads, 1)
            score += min(ratio, 1.0) * 2.0
        if fusion.affected_exons:
            score += 3.5
        if fusion.chrom1 != fusion.chrom2:
            score += 1.5
        if fusion.gene1.upper() in focus_genes:
            score += 2.5
        if fusion.gene2.upper() in focus_genes:
            score += 2.5
        fusion.priority_score = round(score, 3)
        if score >= 16:
            fusion.priority_tier = "Tier1"
        elif score >= 10:
            fusion.priority_tier = "Tier2"
        else:
            fusion.priority_tier = "Tier3"
        if fusion.priority_tier == "Tier1":
            fusion.confidence = "high"
        elif fusion.priority_tier == "Tier2" and fusion.confidence == "low":
            fusion.confidence = "medium"
    return fusions
