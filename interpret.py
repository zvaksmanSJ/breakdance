#!/usr/bin/env python3
"""
interpret.py

Purpose
-------
Create first-pass structural and biological interpretations from
annotated harmonized junctions and clusters.

Version 1 event classes
-----------------------
Junction-driven:
- fusion_candidate
- exon_disruption
- promoter_hijack_candidate
- enhancer_hijack_candidate
- interchromosomal_rearrangement

Cluster-driven:
- complex_rearrangement_cluster
- chromothripsis_like_cluster

Important note
--------------
These are prioritization / review heuristics, not definitive truth labels.
They are intended to help you triage candidates for manual review and
downstream interpretation.
"""

from __future__ import annotations

from schema import InterpretedEvent


def interpret_junctions(sample: str, junctions) -> list[InterpretedEvent]:
    """
    Interpret individual junctions into candidate event types.
    """
    events: list[InterpretedEvent] = []
    event_index = 1

    for j in junctions:
        genes1 = j.annotations.get("genes_side1", [])
        genes2 = j.annotations.get("genes_side2", [])

        exons1 = j.annotations.get("exons_side1", [])
        exons2 = j.annotations.get("exons_side2", [])

        promoters1 = j.annotations.get("promoters_side1", [])
        promoters2 = j.annotations.get("promoters_side2", [])

        enhancers1 = j.annotations.get("enhancers_side1", [])
        enhancers2 = j.annotations.get("enhancers_side2", [])

        # Candidate gene fusion: both sides touch genes and not the exact same set.
        if genes1 and genes2 and set(genes1) != set(genes2):
            events.append(
                InterpretedEvent(
                    event_id=f"EVT_{event_index:06d}",
                    event_type="fusion_candidate",
                    sample=sample,
                    junction_ids=[j.junction_id],
                    genes=sorted(set(genes1 + genes2)),
                    details={
                        "genes_side1": genes1,
                        "genes_side2": genes2,
                        "svtype": j.svtype,
                        "chrom_pair": j.chrom_pair,
                    },
                    confidence=j.confidence,
                )
            )
            event_index += 1

        # Exon disruption: one or both breakends land in exons.
        if exons1 or exons2:
            events.append(
                InterpretedEvent(
                    event_id=f"EVT_{event_index:06d}",
                    event_type="exon_disruption",
                    sample=sample,
                    junction_ids=[j.junction_id],
                    genes=sorted(set(genes1 + genes2)),
                    details={
                        "exons_side1": exons1,
                        "exons_side2": exons2,
                        "svtype": j.svtype,
                    },
                    confidence=j.confidence,
                )
            )
            event_index += 1

        # Promoter hijack candidate: promoter on one side, gene on the other.
        if (promoters1 and genes2) or (promoters2 and genes1):
            events.append(
                InterpretedEvent(
                    event_id=f"EVT_{event_index:06d}",
                    event_type="promoter_hijack_candidate",
                    sample=sample,
                    junction_ids=[j.junction_id],
                    genes=sorted(set(genes1 + genes2)),
                    details={
                        "promoters_side1": promoters1,
                        "promoters_side2": promoters2,
                        "genes_side1": genes1,
                        "genes_side2": genes2,
                    },
                    confidence=j.confidence,
                )
            )
            event_index += 1

        # Enhancer hijack candidate: enhancer on one side, gene on the other.
        if (enhancers1 and genes2) or (enhancers2 and genes1):
            events.append(
                InterpretedEvent(
                    event_id=f"EVT_{event_index:06d}",
                    event_type="enhancer_hijack_candidate",
                    sample=sample,
                    junction_ids=[j.junction_id],
                    genes=sorted(set(genes1 + genes2)),
                    details={
                        "enhancers_side1": enhancers1,
                        "enhancers_side2": enhancers2,
                        "genes_side1": genes1,
                        "genes_side2": genes2,
                    },
                    confidence=j.confidence,
                )
            )
            event_index += 1

        # Generic interchromosomal rearrangement.
        if j.is_interchromosomal:
            events.append(
                InterpretedEvent(
                    event_id=f"EVT_{event_index:06d}",
                    event_type="interchromosomal_rearrangement",
                    sample=sample,
                    junction_ids=[j.junction_id],
                    genes=sorted(set(genes1 + genes2)),
                    details={
                        "chrom_pair": j.chrom_pair,
                        "svtype": j.svtype,
                    },
                    confidence=j.confidence,
                )
            )
            event_index += 1

    return events


def interpret_clusters(sample: str, clusters) -> list[InterpretedEvent]:
    """
    Interpret clusters into larger-scale complex-event candidates.
    """
    events: list[InterpretedEvent] = []
    event_index = 1

    for c in clusters:
        genes = c.annotations.get("genes", [])

        # Simple chromothripsis-like heuristic:
        # many junctions, relatively few chromosomes.
        if c.junction_count >= 8 and len(c.chromosomes) <= 3:
            events.append(
                InterpretedEvent(
                    event_id=f"CLUSTEVT_{event_index:06d}",
                    event_type="chromothripsis_like_cluster",
                    sample=sample,
                    cluster_id=c.cluster_id,
                    junction_ids=c.junction_ids,
                    genes=genes,
                    details={
                        "junction_count": c.junction_count,
                        "chromosomes": c.chromosomes,
                        "interchromosomal_junction_count": c.interchromosomal_junction_count,
                        "notes": c.notes,
                    },
                    confidence="medium",
                )
            )
            event_index += 1

        # Generic complex rearrangement cluster.
        elif c.junction_count >= 5:
            events.append(
                InterpretedEvent(
                    event_id=f"CLUSTEVT_{event_index:06d}",
                    event_type="complex_rearrangement_cluster",
                    sample=sample,
                    cluster_id=c.cluster_id,
                    junction_ids=c.junction_ids,
                    genes=genes,
                    details={
                        "junction_count": c.junction_count,
                        "chromosomes": c.chromosomes,
                        "interchromosomal_junction_count": c.interchromosomal_junction_count,
                        "notes": c.notes,
                    },
                    confidence="medium",
                )
            )
            event_index += 1

    return events
