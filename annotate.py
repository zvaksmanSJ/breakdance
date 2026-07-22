#!/usr/bin/env python3
"""
annotate.py

Add breakpoint and span-based interval annotations to harmonized junctions and clusters.
"""

from __future__ import annotations

from collections import defaultdict


def load_bed4(path: str) -> dict[str, list[tuple[int, int, str]]]:
    intervals = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError(f"BED4 required for annotation file: {path}")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            intervals[chrom].append((start, end, name))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def point_hits(intervals: dict[str, list[tuple[int, int, str]]], chrom: str, pos: int) -> list[str]:
    hits = []
    for start, end, name in intervals.get(chrom, []):
        if start <= pos <= end:
            hits.append(name)
    return sorted(set(hits))


def span_hits(intervals: dict[str, list[tuple[int, int, str]]], chrom: str, start_pos: int, end_pos: int) -> list[str]:
    lo, hi = sorted((start_pos, end_pos))
    hits = []
    for start, end, name in intervals.get(chrom, []):
        if end < lo or start > hi:
            continue
        hits.append(name)
    return sorted(set(hits))


def _strip_feature_names(values: list[str]) -> list[str]:
    cleaned = []
    for value in values:
        cleaned.append(value.split("|")[0])
    return sorted(set(cleaned))


def annotate_junctions(
    junctions,
    genes_bed: str | None = None,
    exons_bed: str | None = None,
    promoters_bed: str | None = None,
    enhancers_bed: str | None = None,
):
    genes = load_bed4(genes_bed) if genes_bed else None
    exons = load_bed4(exons_bed) if exons_bed else None
    promoters = load_bed4(promoters_bed) if promoters_bed else None
    enhancers = load_bed4(enhancers_bed) if enhancers_bed else None

    for j in junctions:
        ann = j.annotations

        if genes:
            ann["genes_side1"] = point_hits(genes, j.chrom1, j.pos1)
            ann["genes_side2"] = point_hits(genes, j.chrom2, j.pos2)
            if not j.is_interchromosomal and j.svtype in {"DEL", "DUP", "INV"}:
                ann["genes_span"] = span_hits(genes, j.chrom1, j.pos1, j.pos2)
            else:
                ann["genes_span"] = []

        if exons:
            ann["exons_side1"] = point_hits(exons, j.chrom1, j.pos1)
            ann["exons_side2"] = point_hits(exons, j.chrom2, j.pos2)
            if not j.is_interchromosomal and j.svtype in {"DEL", "DUP", "INV"}:
                ann["exons_span"] = span_hits(exons, j.chrom1, j.pos1, j.pos2)
            else:
                ann["exons_span"] = []

        if promoters:
            ann["promoters_side1"] = point_hits(promoters, j.chrom1, j.pos1)
            ann["promoters_side2"] = point_hits(promoters, j.chrom2, j.pos2)
            if not j.is_interchromosomal and j.svtype in {"DEL", "DUP", "INV"}:
                ann["promoters_span"] = span_hits(promoters, j.chrom1, j.pos1, j.pos2)
            else:
                ann["promoters_span"] = []

        if enhancers:
            ann["enhancers_side1"] = point_hits(enhancers, j.chrom1, j.pos1)
            ann["enhancers_side2"] = point_hits(enhancers, j.chrom2, j.pos2)
            if not j.is_interchromosomal and j.svtype in {"DEL", "DUP", "INV"}:
                ann["enhancers_span"] = span_hits(enhancers, j.chrom1, j.pos1, j.pos2)
            else:
                ann["enhancers_span"] = []

        gene_names = _strip_feature_names(
            ann.get("genes_side1", []) + ann.get("genes_side2", []) + ann.get("genes_span", [])
        )
        exon_names = sorted(set(ann.get("exons_side1", []) + ann.get("exons_side2", []) + ann.get("exons_span", [])))
        promoter_names = sorted(set(ann.get("promoters_side1", []) + ann.get("promoters_side2", []) + ann.get("promoters_span", [])))
        enhancer_names = sorted(set(ann.get("enhancers_side1", []) + ann.get("enhancers_side2", []) + ann.get("enhancers_span", [])))

        j.affected_genes = gene_names
        j.affected_exons = exon_names
        j.affected_promoters = promoter_names
        j.affected_enhancers = enhancer_names
        ann["affected_genes"] = gene_names
        ann["affected_exons"] = exon_names
        ann["affected_promoters"] = promoter_names
        ann["affected_enhancers"] = enhancer_names

    return junctions


def annotate_clusters(clusters, junction_by_id):
    for c in clusters:
        genes = set()
        exons = set()
        promoters = set()
        enhancers = set()
        spans = {}

        for jid in c.junction_ids:
            j = junction_by_id[jid]
            ann = j.annotations
            genes.update(ann.get("affected_genes", []))
            exons.update(ann.get("affected_exons", []))
            promoters.update(ann.get("affected_promoters", []))
            enhancers.update(ann.get("affected_enhancers", []))
            for chrom, pos in ((j.chrom1, j.pos1), (j.chrom2, j.pos2)):
                if chrom not in spans:
                    spans[chrom] = [pos, pos]
                else:
                    spans[chrom][0] = min(spans[chrom][0], pos)
                    spans[chrom][1] = max(spans[chrom][1], pos)

        c.annotations["genes"] = sorted(genes)
        c.annotations["exons"] = sorted(exons)
        c.annotations["promoters"] = sorted(promoters)
        c.annotations["enhancers"] = sorted(enhancers)
        c.annotations["event_span_by_chrom"] = spans
        c.affected_genes = sorted(genes)
        c.affected_exons = sorted(exons)
        c.affected_promoters = sorted(promoters)
        c.affected_enhancers = sorted(enhancers)
        c.event_span_by_chrom = spans

    return clusters
