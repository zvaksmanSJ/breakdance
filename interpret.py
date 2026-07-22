#!/usr/bin/env python3
"""
interpret.py

Create first-pass structural and biological interpretations from
annotated harmonized junctions and clusters, including dedicated fusion events.
"""

from __future__ import annotations

from schema import FusionEvent, InterpretedEvent


def _location_summary(chrom1: str, pos1: int, chrom2: str, pos2: int) -> str:
    if chrom1 == chrom2:
        lo, hi = sorted((pos1, pos2))
        return f"{chrom1}:{lo}-{hi}"
    return f"{chrom1}:{pos1} <-> {chrom2}:{pos2}"


def _confidence_from_callers(caller_count: int, supporting_reads: int | None) -> str:
    if caller_count >= 3:
        return "high"
    if caller_count == 2:
        return "medium"
    if supporting_reads is not None and supporting_reads >= 8:
        return "medium"
    return "low"


def _fusion_candidates(sample: str, junctions) -> tuple[list[FusionEvent], list[InterpretedEvent]]:
    fusion_events: list[FusionEvent] = []
    interpreted: list[InterpretedEvent] = []
    event_index = 1

    for j in junctions:
        genes1 = j.annotations.get("genes_side1", [])
        genes2 = j.annotations.get("genes_side2", [])
        if not genes1 or not genes2:
            continue

        left_genes = sorted({g.split("|")[0] for g in genes1})
        right_genes = sorted({g.split("|")[0] for g in genes2})
        if set(left_genes) == set(right_genes):
            continue

        gene1 = left_genes[0]
        gene2 = right_genes[0]
        label = f"{gene1}--{gene2}"
        confidence = _confidence_from_callers(j.caller_count, j.supporting_reads)
        location_summary = _location_summary(j.chrom1, j.pos1, j.chrom2, j.pos2)
        fusion_id = f"FUSION_{event_index:06d}"

        fusion = FusionEvent(
            fusion_id=fusion_id,
            sample=sample,
            gene1=gene1,
            gene2=gene2,
            chrom1=j.chrom1,
            pos1=j.pos1,
            chrom2=j.chrom2,
            pos2=j.pos2,
            orientation=f"{j.strand1 or '?'}{j.strand2 or '?'}",
            svtype=j.svtype,
            junction_ids=[j.junction_id],
            supporting_callers=list(j.callers),
            support_by_caller=dict(j.support_by_caller),
            caller_count=j.caller_count,
            supporting_reads=j.supporting_reads,
            covering_reads=j.covering_reads,
            non_supporting_reads=j.non_supporting_reads,
            affected_exons=list(j.affected_exons),
            affected_genes=list(j.affected_genes),
            location_summary=location_summary,
            confidence=confidence,
            details={
                "genes_side1": left_genes,
                "genes_side2": right_genes,
                "chrom_pair": j.chrom_pair,
                "interchromosomal": j.is_interchromosomal,
                "event_region1": j.event_region1,
                "event_region2": j.event_region2,
            },
        )
        fusion_events.append(fusion)

        interpreted.append(
            InterpretedEvent(
                event_id=fusion_id,
                event_type="fusion_event",
                sample=sample,
                event_label=label,
                junction_ids=[j.junction_id],
                location_summary=location_summary,
                genes=sorted(set(left_genes + right_genes)),
                affected_genes=list(j.affected_genes),
                supporting_callers=list(j.callers),
                support_by_caller=dict(j.support_by_caller),
                supporting_reads=j.supporting_reads,
                covering_reads=j.covering_reads,
                non_supporting_reads=j.non_supporting_reads,
                details={
                    "fusion_id": fusion_id,
                    "gene1": gene1,
                    "gene2": gene2,
                    "orientation": fusion.orientation,
                    "affected_exons": list(j.affected_exons),
                },
                confidence=confidence,
            )
        )
        event_index += 1

    return fusion_events, interpreted


def interpret_junctions(sample: str, junctions) -> tuple[list[FusionEvent], list[InterpretedEvent]]:
    fusion_events, events = _fusion_candidates(sample, junctions)
    event_index = 1

    for j in junctions:
        genes = list(j.affected_genes)
        location_summary = _location_summary(j.chrom1, j.pos1, j.chrom2, j.pos2)
        confidence = _confidence_from_callers(j.caller_count, j.supporting_reads)

        if j.affected_exons:
            events.append(
                InterpretedEvent(
                    event_id=f"EXON_{event_index:06d}",
                    event_type="exon_interfering_sv",
                    sample=sample,
                    event_label="exon_interference",
                    junction_ids=[j.junction_id],
                    location_summary=location_summary,
                    genes=genes,
                    affected_genes=genes,
                    supporting_callers=list(j.callers),
                    support_by_caller=dict(j.support_by_caller),
                    supporting_reads=j.supporting_reads,
                    covering_reads=j.covering_reads,
                    non_supporting_reads=j.non_supporting_reads,
                    details={
                        "svtype": j.svtype,
                        "affected_exons": list(j.affected_exons),
                        "event_region1": j.event_region1,
                        "event_region2": j.event_region2,
                    },
                    confidence=confidence,
                )
            )
            event_index += 1

        if (j.annotations.get("promoters_side1") and j.annotations.get("genes_side2")) or (
            j.annotations.get("promoters_side2") and j.annotations.get("genes_side1")
        ):
            events.append(
                InterpretedEvent(
                    event_id=f"PROM_{event_index:06d}",
                    event_type="promoter_hijack_candidate",
                    sample=sample,
                    event_label="promoter_hijack",
                    junction_ids=[j.junction_id],
                    location_summary=location_summary,
                    genes=genes,
                    affected_genes=genes,
                    supporting_callers=list(j.callers),
                    support_by_caller=dict(j.support_by_caller),
                    supporting_reads=j.supporting_reads,
                    details={
                        "promoters": list(j.affected_promoters),
                        "genes_side1": j.annotations.get("genes_side1", []),
                        "genes_side2": j.annotations.get("genes_side2", []),
                    },
                    confidence=confidence,
                )
            )
            event_index += 1

        if (j.annotations.get("enhancers_side1") and j.annotations.get("genes_side2")) or (
            j.annotations.get("enhancers_side2") and j.annotations.get("genes_side1")
        ):
            events.append(
                InterpretedEvent(
                    event_id=f"ENH_{event_index:06d}",
                    event_type="enhancer_hijack_candidate",
                    sample=sample,
                    event_label="enhancer_hijack",
                    junction_ids=[j.junction_id],
                    location_summary=location_summary,
                    genes=genes,
                    affected_genes=genes,
                    supporting_callers=list(j.callers),
                    support_by_caller=dict(j.support_by_caller),
                    supporting_reads=j.supporting_reads,
                    details={
                        "enhancers": list(j.affected_enhancers),
                        "genes_side1": j.annotations.get("genes_side1", []),
                        "genes_side2": j.annotations.get("genes_side2", []),
                    },
                    confidence=confidence,
                )
            )
            event_index += 1

        if j.is_interchromosomal:
            events.append(
                InterpretedEvent(
                    event_id=f"CTX_{event_index:06d}",
                    event_type="interchromosomal_rearrangement",
                    sample=sample,
                    event_label="interchromosomal_event",
                    junction_ids=[j.junction_id],
                    location_summary=location_summary,
                    genes=genes,
                    affected_genes=genes,
                    supporting_callers=list(j.callers),
                    support_by_caller=dict(j.support_by_caller),
                    supporting_reads=j.supporting_reads,
                    details={
                        "chrom_pair": j.chrom_pair,
                        "svtype": j.svtype,
                    },
                    confidence=confidence,
                )
            )
            event_index += 1
        elif j.svtype == "DEL" and len(j.affected_genes) >= 1:
            events.append(
                InterpretedEvent(
                    event_id=f"DEL_{event_index:06d}",
                    event_type="large_deletion",
                    sample=sample,
                    event_label="large_deletion",
                    junction_ids=[j.junction_id],
                    location_summary=location_summary,
                    genes=genes,
                    affected_genes=genes,
                    supporting_callers=list(j.callers),
                    support_by_caller=dict(j.support_by_caller),
                    supporting_reads=j.supporting_reads,
                    details={
                        "svtype": j.svtype,
                        "affected_genes": list(j.affected_genes),
                        "event_region1": j.event_region1,
                        "event_region2": j.event_region2,
                    },
                    confidence=confidence,
                )
            )
            event_index += 1

    return fusion_events, events


def interpret_clusters(sample: str, clusters) -> list[InterpretedEvent]:
    events: list[InterpretedEvent] = []
    event_index = 1

    for c in clusters:
        genes = list(c.affected_genes)
        location_summary = "; ".join(
            f"{chrom}:{span[0]}-{span[1]}" for chrom, span in sorted(c.event_span_by_chrom.items())
        )

        if c.junction_count >= 8 and len(c.chromosomes) <= 3:
            events.append(
                InterpretedEvent(
                    event_id=f"CLUST_{event_index:06d}",
                    event_type="chromothripsis_like_cluster",
                    sample=sample,
                    event_label="chromothripsis_like_cluster",
                    cluster_id=c.cluster_id,
                    junction_ids=c.junction_ids,
                    location_summary=location_summary,
                    genes=genes,
                    affected_genes=genes,
                    supporting_callers=list(c.caller_set),
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
        elif c.junction_count >= 5:
            events.append(
                InterpretedEvent(
                    event_id=f"CLUST_{event_index:06d}",
                    event_type="complex_rearrangement_cluster",
                    sample=sample,
                    event_label="complex_rearrangement_cluster",
                    cluster_id=c.cluster_id,
                    junction_ids=c.junction_ids,
                    location_summary=location_summary,
                    genes=genes,
                    affected_genes=genes,
                    supporting_callers=list(c.caller_set),
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
