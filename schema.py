#!/usr/bin/env python3
"""
schema.py

Core dataclasses for Breakdance multi-caller SV interpretation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Optional


@dataclass
class RawSVCall:
    caller: str
    record_id: str
    svtype: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    strand1: Optional[str] = None
    strand2: Optional[str] = None
    svlen: Optional[int] = None
    support: Optional[int] = None
    qual: Optional[float] = None
    filt: str = "PASS"
    source_path: Optional[str] = None
    sample_name: Optional[str] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    filter_status: Optional[str] = None
    event_length: Optional[int] = None
    annotations: dict[str, Any] = field(default_factory=dict)
    info: dict[str, Any] = field(default_factory=dict)
    adaptive_local_density_hint: Optional[float] = None
    adaptive_hotspot_hint: Optional[bool] = None


@dataclass
class HarmonizedJunction:
    junction_id: str
    svtype: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    strand1: Optional[str] = None
    strand2: Optional[str] = None
    callers: list[str] = field(default_factory=list)
    raw_call_ids: list[str] = field(default_factory=list)
    support_by_caller: dict[str, Optional[int]] = field(default_factory=dict)
    qual_by_caller: dict[str, Optional[float]] = field(default_factory=dict)
    filters_by_caller: dict[str, str] = field(default_factory=dict)
    source_files: list[str] = field(default_factory=list)
    caller_count: int = 0
    median_support: Optional[float] = None
    confidence: str = "low"
    chrom_pair: Optional[str] = None
    is_interchromosomal: bool = False
    adaptive_hotspot_cluster: bool = False
    adaptive_density_score: Optional[float] = None
    consensus_event_type: Optional[str] = None
    event_region1: Optional[str] = None
    event_region2: Optional[str] = None
    affected_genes: list[str] = field(default_factory=list)
    affected_exons: list[str] = field(default_factory=list)
    affected_promoters: list[str] = field(default_factory=list)
    affected_enhancers: list[str] = field(default_factory=list)
    supporting_reads: Optional[int] = None
    covering_reads: Optional[int] = None
    non_supporting_reads: Optional[int] = None
    annotations: dict[str, Any] = field(default_factory=dict)


@dataclass
class SVCluster:
    cluster_id: str
    junction_ids: list[str] = field(default_factory=list)
    chromosomes: list[str] = field(default_factory=list)
    junction_count: int = 0
    interchromosomal_junction_count: int = 0
    caller_set: list[str] = field(default_factory=list)
    max_adaptive_density_score: Optional[float] = None
    any_adaptive_hotspot: bool = False
    min_position_by_chrom: dict[str, int] = field(default_factory=dict)
    max_position_by_chrom: dict[str, int] = field(default_factory=dict)
    event_span_by_chrom: dict[str, list[int]] = field(default_factory=dict)
    affected_genes: list[str] = field(default_factory=list)
    affected_exons: list[str] = field(default_factory=list)
    affected_promoters: list[str] = field(default_factory=list)
    affected_enhancers: list[str] = field(default_factory=list)
    figure_path: Optional[str] = None
    priority_score: float = 0.0
    priority_tier: str = "Tier3"
    notes: dict[str, Any] = field(default_factory=dict)
    annotations: dict[str, Any] = field(default_factory=dict)


@dataclass
class FusionEvent:
    fusion_id: str
    sample: str
    gene1: str
    gene2: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    orientation: Optional[str] = None
    svtype: str = "BND"
    junction_ids: list[str] = field(default_factory=list)
    supporting_callers: list[str] = field(default_factory=list)
    support_by_caller: dict[str, Optional[int]] = field(default_factory=dict)
    caller_count: int = 0
    supporting_reads: Optional[int] = None
    covering_reads: Optional[int] = None
    non_supporting_reads: Optional[int] = None
    affected_exons: list[str] = field(default_factory=list)
    affected_genes: list[str] = field(default_factory=list)
    location_summary: str = ""
    figure_path: Optional[str] = None
    confidence: str = "low"
    priority_score: float = 0.0
    priority_tier: str = "Tier3"
    details: dict[str, Any] = field(default_factory=dict)

    @property
    def fusion_label(self) -> str:
        return f"{self.gene1}--{self.gene2}"


@dataclass
class InterpretedEvent:
    event_id: str
    event_type: str
    sample: str
    event_label: Optional[str] = None
    junction_ids: list[str] = field(default_factory=list)
    cluster_id: Optional[str] = None
    location_summary: str = ""
    genes: list[str] = field(default_factory=list)
    affected_genes: list[str] = field(default_factory=list)
    supporting_callers: list[str] = field(default_factory=list)
    support_by_caller: dict[str, Optional[int]] = field(default_factory=dict)
    supporting_reads: Optional[int] = None
    covering_reads: Optional[int] = None
    non_supporting_reads: Optional[int] = None
    figure_paths: list[str] = field(default_factory=list)
    details: dict[str, Any] = field(default_factory=dict)
    confidence: str = "low"
    priority_score: float = 0.0
    priority_tier: str = "Tier3"


@dataclass
class BreakpointGraph:
    sample: str
    nodes: list[dict[str, Any]] = field(default_factory=list)
    edges: list[dict[str, Any]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


def dataclass_to_dict(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        return asdict(value)
    return value
