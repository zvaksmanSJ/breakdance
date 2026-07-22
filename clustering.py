#!/usr/bin/env python3
"""
clustering.py

Purpose
-------
Cluster harmonized breakpoint junctions into candidate event groups.

Why this module exists
----------------------
A sample may contain many harmonized junctions, but not all of them belong
to independent events. Some are part of the same local rearrangement,
translocation chain, inversion cluster, or chromothripsis-like locus.

This module groups harmonized junctions into clusters using a simple,
transparent proximity rule.

Version 1 clustering rule
-------------------------
Two junctions are connected if any breakpoint from junction A is within
`distance_bp` of any breakpoint from junction B on the same chromosome.

This is intentionally simple and review-friendly.

What this module does
---------------------
- groups harmonized junctions into connected components
- tracks chromosomes represented in each cluster
- tracks interchromosomal content
- propagates adaptive-density hints from harmonized junctions

What this module does not do
----------------------------
- no graph-theoretic segment reconstruction
- no orientation-aware chaining
- no copy-number integration
"""

from __future__ import annotations

from collections import defaultdict

from schema import HarmonizedJunction, SVCluster


class UnionFind:
    """
    Simple disjoint-set / union-find structure.

    Used to compute connected components among junctions.
    """

    def __init__(self, items: list[str]):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x: str) -> str:
        """
        Return the root representative for x.
        """
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        """
        Merge the sets containing a and b.
        """
        ra = self.find(a)
        rb = self.find(b)

        if ra == rb:
            return

        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def breakpoint_pairs_for_junction(j: HarmonizedJunction) -> list[tuple[str, int]]:
    """
    Return the two breakpoint coordinates for a harmonized junction.
    """
    return [
        (j.chrom1, j.pos1),
        (j.chrom2, j.pos2),
    ]


def junctions_are_clusterable(
    a: HarmonizedJunction,
    b: HarmonizedJunction,
    distance_bp: int,
) -> bool:
    """
    Determine whether two harmonized junctions should be placed in the same cluster.

    Rule
    ----
    If any breakpoint from A lies within `distance_bp` of any breakpoint
    from B on the same chromosome, connect them.
    """
    a_points = breakpoint_pairs_for_junction(a)
    b_points = breakpoint_pairs_for_junction(b)

    for a_chrom, a_pos in a_points:
        for b_chrom, b_pos in b_points:
            if a_chrom != b_chrom:
                continue
            if abs(a_pos - b_pos) <= distance_bp:
                return True

    return False


def cluster_junctions(
    junctions: list[HarmonizedJunction],
    distance_bp: int = 10_000,
) -> list[SVCluster]:
    """
    Cluster harmonized junctions into connected components.

    Parameters
    ----------
    junctions
        Harmonized junctions for one sample.

    distance_bp
        Maximum same-chromosome breakpoint distance for linking two junctions.

    Returns
    -------
    List of SVCluster objects.

    Notes
    -----
    This Version 1 implementation uses O(n^2) pairwise comparison, which is
    acceptable for prototype-scale per-sample junction sets.
    """
    if not junctions:
        return []

    by_id = {j.junction_id: j for j in junctions}
    ids = [j.junction_id for j in junctions]

    uf = UnionFind(ids)

    for i in range(len(junctions)):
        for k in range(i + 1, len(junctions)):
            a = junctions[i]
            b = junctions[k]

            if junctions_are_clusterable(a, b, distance_bp=distance_bp):
                uf.union(a.junction_id, b.junction_id)

    grouped_ids = defaultdict(list)
    for jid in ids:
        grouped_ids[uf.find(jid)].append(jid)

    clusters: list[SVCluster] = []

    for idx, member_ids in enumerate(grouped_ids.values(), start=1):
        members = [by_id[jid] for jid in member_ids]

        chroms = set()
        callers = set()
        interchrom_count = 0
        adaptive_scores = []
        any_hotspot = False

        min_pos = {}
        max_pos = {}

        for j in members:
            chroms.add(j.chrom1)
            chroms.add(j.chrom2)
            callers.update(j.callers)

            if j.is_interchromosomal:
                interchrom_count += 1

            if j.adaptive_density_score is not None:
                adaptive_scores.append(j.adaptive_density_score)

            if j.adaptive_hotspot_cluster:
                any_hotspot = True

            for chrom, pos in ((j.chrom1, j.pos1), (j.chrom2, j.pos2)):
                if chrom not in min_pos or pos < min_pos[chrom]:
                    min_pos[chrom] = pos
                if chrom not in max_pos or pos > max_pos[chrom]:
                    max_pos[chrom] = pos

        notes = {}

        if len(member_ids) >= 5:
            notes["complexity_hint"] = "multi_junction_cluster"

        if interchrom_count > 0:
            notes["interchromosomal_hint"] = "contains_interchromosomal_junctions"

        if any_hotspot:
            notes["adaptive_density_hint"] = "dense_breakpoint_region"

        clusters.append(
            SVCluster(
                cluster_id=f"CLUSTER_{idx:06d}",
                junction_ids=sorted(member_ids),
                chromosomes=sorted(chroms),
                junction_count=len(member_ids),
                interchromosomal_junction_count=interchrom_count,
                caller_set=sorted(callers),
                max_adaptive_density_score=max(adaptive_scores) if adaptive_scores else None,
                any_adaptive_hotspot=any_hotspot,
                min_position_by_chrom=min_pos,
                max_position_by_chrom=max_pos,
                notes=notes,
            )
        )

    # Order larger / potentially more interesting clusters first.
    clusters.sort(
        key=lambda c: (
            -c.junction_count,
            -c.interchromosomal_junction_count,
            c.cluster_id,
        )
    )

    return clusters
