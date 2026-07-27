#!/usr/bin/env python3
"""
coverage.py

Optional read-evidence helpers for Breakdance.

This module is intentionally lightweight and degrades gracefully when BAM
inputs or pysam are unavailable. Current behavior focuses on surfacing basic
coverage-derived context rather than full breakpoint-spanning read validation.
"""

from __future__ import annotations

from pathlib import Path

try:
    import pysam  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    pysam = None


DEFAULT_WINDOW_BP = 250


def _safe_mean(values: list[int]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def _estimate_depth_near_breakpoint(bam_path: Path, chrom: str, pos: int, window_bp: int) -> int | None:
    if pysam is None:
        return None
    if not bam_path.exists():
        return None

    start = max(0, pos - window_bp)
    end = max(start + 1, pos + window_bp)

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            depths: list[int] = []
            for column in bam.pileup(chrom, start, end, truncate=True, stepper="nofilter"):
                if start <= column.reference_pos < end:
                    depths.append(column.nsegments)
    except Exception:
        return None

    mean_depth = _safe_mean(depths)
    if mean_depth is None:
        return None
    return int(round(mean_depth))


def attach_read_evidence(junctions, bam_path: str | Path | None, window_bp: int = DEFAULT_WINDOW_BP):
    """
    Populate approximate coverage-related fields on harmonized junctions.

    Notes
    -----
    - `supporting_reads` is preserved if already set from caller support.
    - `covering_reads` is approximated from local depth around both breakpoints.
    - `non_supporting_reads` is approximated as max(covering - supporting, 0).
    - If BAM support is unavailable, junctions are returned unchanged.
    """
    if bam_path is None:
        return junctions

    bam_path = Path(bam_path)
    if pysam is None or not bam_path.exists():
        return junctions

    for j in junctions:
        depth1 = _estimate_depth_near_breakpoint(bam_path, j.chrom1, j.pos1, window_bp)
        depth2 = _estimate_depth_near_breakpoint(bam_path, j.chrom2, j.pos2, window_bp)
        depth_values = [d for d in (depth1, depth2) if d is not None]

        if depth_values:
            j.covering_reads = int(round(sum(depth_values) / len(depth_values)))
            if j.supporting_reads is not None:
                j.non_supporting_reads = max(j.covering_reads - j.supporting_reads, 0)

        j.annotations.setdefault("read_evidence", {})
        j.annotations["read_evidence"].update(
            {
                "bam_path": str(bam_path),
                "window_bp": window_bp,
                "depth_near_side1": depth1,
                "depth_near_side2": depth2,
                "coverage_estimated": bool(depth_values),
            }
        )

    return junctions


def propagate_read_evidence_to_events(events, junction_by_id):
    for e in events:
        supports = []
        coverings = []
        non_supportings = []
        for jid in e.junction_ids:
            j = junction_by_id.get(jid)
            if not j:
                continue
            if j.supporting_reads is not None:
                supports.append(j.supporting_reads)
            if j.covering_reads is not None:
                coverings.append(j.covering_reads)
            if j.non_supporting_reads is not None:
                non_supportings.append(j.non_supporting_reads)

        if supports:
            e.supporting_reads = max(supports)
        if coverings:
            e.covering_reads = max(coverings)
        if non_supportings:
            e.non_supporting_reads = min(non_supportings)
    return events


def propagate_read_evidence_to_fusions(fusions, junction_by_id):
    for fusion in fusions:
        supports = []
        coverings = []
        non_supportings = []
        for jid in fusion.junction_ids:
            j = junction_by_id.get(jid)
            if not j:
                continue
            if j.supporting_reads is not None:
                supports.append(j.supporting_reads)
            if j.covering_reads is not None:
                coverings.append(j.covering_reads)
            if j.non_supporting_reads is not None:
                non_supportings.append(j.non_supporting_reads)

        if supports:
            fusion.supporting_reads = max(supports)
        if coverings:
            fusion.covering_reads = max(coverings)
        if non_supportings:
            fusion.non_supporting_reads = min(non_supportings)
    return fusions
