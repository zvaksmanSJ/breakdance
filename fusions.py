#!/usr/bin/env python3
"""
fusions.py

Helpers for deduplicating and normalizing fusion event outputs.
"""

from __future__ import annotations

from collections import defaultdict


FUSION_BIN_BP = 1_000


def _bin_pos(pos: int) -> int:
    return int(pos) // FUSION_BIN_BP


def _fusion_key(fusion):
    genes = tuple(sorted([fusion.gene1, fusion.gene2]))
    chroms = tuple(sorted([fusion.chrom1, fusion.chrom2]))
    pos_bins = tuple(sorted([_bin_pos(fusion.pos1), _bin_pos(fusion.pos2)]))
    orientation = fusion.orientation or "??"
    return genes, chroms, pos_bins, orientation


def consolidate_fusions(fusions):
    grouped = defaultdict(list)
    for fusion in fusions:
        grouped[_fusion_key(fusion)].append(fusion)

    consolidated = []
    for group in grouped.values():
        best = sorted(
            group,
            key=lambda f: (
                -float(f.priority_score),
                -int(f.caller_count),
                f.fusion_id,
            ),
        )[0]

        all_junction_ids = sorted({jid for f in group for jid in f.junction_ids})
        all_callers = sorted({caller for f in group for caller in f.supporting_callers})
        support_by_caller = {}
        for f in group:
            support_by_caller.update(f.support_by_caller)

        best.junction_ids = all_junction_ids
        best.supporting_callers = all_callers
        best.support_by_caller = support_by_caller
        best.affected_exons = sorted({x for f in group for x in f.affected_exons})
        best.affected_genes = sorted({x for f in group for x in f.affected_genes})
        best.details = dict(best.details)
        best.details["consolidated_fusion_count"] = len(group)
        best.details["fusion_bin_bp"] = FUSION_BIN_BP
        consolidated.append(best)

    consolidated.sort(key=lambda f: (-float(f.priority_score), f.fusion_label, f.fusion_id))
    return consolidated
