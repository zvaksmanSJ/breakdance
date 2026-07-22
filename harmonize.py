#!/usr/bin/env python3
from __future__ import annotations
from collections import defaultdict
from statistics import median
from schema import HarmonizedJunction, RawSVCall

def canonicalize_call(call: RawSVCall) -> RawSVCall:
    left = (call.chrom1, call.pos1)
    right = (call.chrom2, call.pos2)
    if left <= right:
        return call
    return RawSVCall(
        caller=call.caller,
        record_id=call.record_id,
        svtype=call.svtype,
        chrom1=call.chrom2,
        pos1=call.pos2,
        chrom2=call.chrom1,
        pos2=call.pos1,
        strand1=call.strand2,
        strand2=call.strand1,
        svlen=call.svlen,
        support=call.support,
        qual=call.qual,
        filt=call.filt,
        source_path=call.source_path,
        sample_name=call.sample_name,
        info=call.info,
        adaptive_local_density_hint=call.adaptive_local_density_hint,
        adaptive_hotspot_hint=call.adaptive_hotspot_hint,
    )

def compatible_svtypes(a: str, b: str) -> bool:
    if a == b:
        return True
    if {a, b} <= {"BND"}:
        return True
    return False

def calls_match(a: RawSVCall, b: RawSVCall, tolerance_bp: int) -> bool:
    if not compatible_svtypes(a.svtype, b.svtype):
        return False
    if a.chrom1 != b.chrom1 or a.chrom2 != b.chrom2:
        return False
    if abs(a.pos1 - b.pos1) > tolerance_bp:
        return False
    if abs(a.pos2 - b.pos2) > tolerance_bp:
        return False
    return True

def annotate_density_hints(raw_calls: list[RawSVCall], window_bp: int = 1000000) -> list[RawSVCall]:
    by_chrom = defaultdict(list)
    for idx, call in enumerate(raw_calls):
        by_chrom[call.chrom1].append((call.pos1, idx))
        by_chrom[call.chrom2].append((call.pos2, idx))

    density_scores = [0.0] * len(raw_calls)

    for chrom, entries in by_chrom.items():
        entries.sort()
        positions = [x[0] for x in entries]
        left = 0
        for right in range(len(entries)):
            while positions[right] - positions[left] > window_bp:
                left += 1
            window_size = right - left + 1
            for j in range(left, right + 1):
                idx = entries[j][1]
                density_scores[idx] = max(density_scores[idx], float(window_size))

    unique_scores = sorted(set(density_scores)) if density_scores else [0.0]
    hotspot_threshold = unique_scores[max(0, int(0.95 * (len(unique_scores) - 1)))] if unique_scores else 0.0

    for i, call in enumerate(raw_calls):
        call.adaptive_local_density_hint = density_scores[i]
        call.adaptive_hotspot_hint = density_scores[i] >= hotspot_threshold and density_scores[i] > 0

    return raw_calls

def harmonize_calls(raw_calls: list[RawSVCall], tolerance_bp: int = 100) -> list[HarmonizedJunction]:
    canon_calls = [canonicalize_call(c) for c in raw_calls]
    canon_calls.sort(key=lambda x: (x.chrom1, x.pos1, x.chrom2, x.pos2, x.svtype, x.caller))

    groups = []
    for call in canon_calls:
        placed = False
        for group in groups:
            if calls_match(call, group[0], tolerance_bp=tolerance_bp):
                group.append(call)
                placed = True
                break
        if not placed:
            groups.append([call])

    harmonized = []
    for i, group in enumerate(groups, start=1):
        rep = group[0]
        callers = sorted({c.caller for c in group})
        supports = [c.support for c in group if c.support is not None]
        support_by_caller = {c.caller: c.support for c in group}
        qual_by_caller = {c.caller: c.qual for c in group}
        filters_by_caller = {c.caller: c.filt for c in group}
        raw_call_ids = [f"{c.caller}:{c.record_id}" for c in group]
        source_files = sorted({str(c.source_path) for c in group})
        median_support = median(supports) if supports else None
        caller_count = len(callers)

        confidence = "low"
        if caller_count >= 3:
            confidence = "high"
        elif caller_count == 2:
            confidence = "medium"
        elif caller_count == 1 and median_support is not None and median_support >= 5:
            confidence = "medium"

        density_scores = [c.adaptive_local_density_hint for c in group if c.adaptive_local_density_hint is not None]
        density_score = max(density_scores) if density_scores else None
        hotspot_cluster = any(bool(c.adaptive_hotspot_hint) for c in group)

        harmonized.append(
            HarmonizedJunction(
                junction_id=f"JUNC_{i:06d}",
                svtype=rep.svtype,
                chrom1=rep.chrom1,
                pos1=round(sum(c.pos1 for c in group) / len(group)),
                chrom2=rep.chrom2,
                pos2=round(sum(c.pos2 for c in group) / len(group)),
                callers=callers,
                raw_call_ids=raw_call_ids,
                support_by_caller=support_by_caller,
                qual_by_caller=qual_by_caller,
                filters_by_caller=filters_by_caller,
                source_files=source_files,
                caller_count=caller_count,
                median_support=median_support,
                confidence=confidence,
                chrom_pair=f"{rep.chrom1}|{rep.chrom2}",
                is_interchromosomal=(rep.chrom1 != rep.chrom2),
                adaptive_hotspot_cluster=hotspot_cluster,
                adaptive_density_score=density_score,
            )
        )
    return harmonized
