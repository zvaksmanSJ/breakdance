#!/usr/bin/env python3
"""
clustering.py

Purpose
-------
Cluster harmonized breakpoint junctions into candidate event groups.
"""

from __future__ import annotations

from collections import defaultdict

from schema import HarmonizedJunction, SVCluster


class UnionFind:
    def __init__(self, items: list[str]):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x: str) -> str:
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
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
    return [(j.chrom1, j.pos1), (j.chrom2, j.pos2)]


def junctions_are_clusterable(a: HarmonizedJunction, b: HarmonizedJunction, distance_bp: int) -> bool:
    a_points = breakpoint_pairs_for_junction(a)
    b_points = breakpoint_pairs_for_junction(b)
    for a_chrom, a_pos in a_points:
        for b_chrom, b_pos in b_points:
            if a_chrom != b_chrom:
                continue
            if abs(a_pos - b_pos) <= distance_bp:
                return True
    return False


def cluster_junctions(junctions: list[HarmonizedJunction], distance_bp: int = 10_000) -> list[SVCluster]:
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
        affected_genes = set()
        affected_exons = set()
        affected_promoters = set()
        affected_enhancers = set()

        for j in members:
            chroms.add(j.chrom1)
            chroms.add(j.chrom2)
            callers.update(j.callers)
            affected_genes.update(j.affected_genes)
            affected_exons.update(j.affected_exons)
            affected_promoters.update(j.affected_promoters)
            affected_enhancers.update(j.affected_enhancers)

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

        event_span_by_chrom = {
            chrom: [min_pos[chrom], max_pos[chrom]]
            for chrom in sorted(min_pos)
        }

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
                event_span_by_chrom=event_span_by_chrom,
                affected_genes=sorted(affected_genes),
                affected_exons=sorted(affected_exons),
                affected_promoters=sorted(affected_promoters),
                affected_enhancers=sorted(affected_enhancers),
                notes=notes,
            )
        )

    clusters.sort(key=lambda c: (-c.junction_count, -c.interchromosomal_junction_count, c.cluster_id))
    return clusters
