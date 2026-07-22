#!/usr/bin/env python3
"""
annotate.py

Purpose
-------
Annotate harmonized junctions and clusters with genomic interval features.

Supported interval sets
-----------------------
- genes
- exons
- promoters
- enhancers

Input format
------------
Each interval file is expected to be BED4:
    chrom    start    end    name

Example
-------
chr8    127735433    127742951    MYC
chr14   105586337    105864198    IGH

Design notes
------------
This starter version uses a simple chromosome-keyed list scan.
That keeps the code very readable and easy to debug.

For very large annotation files, this can later be replaced with interval trees.
"""

from __future__ import annotations

from collections import defaultdict


def load_bed4(path: str) -> dict[str, list[tuple[int, int, str]]]:
    """
    Load a BED4 file into a chromosome-indexed interval dictionary.

    Returns
    -------
    dict:
        chrom -> list of (start, end, name)
    """
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

    return intervals


def point_hits(
    intervals: dict[str, list[tuple[int, int, str]]],
    chrom: str,
    pos: int,
) -> list[str]:
    """
    Return all interval names overlapping a single breakpoint coordinate.

    Notes
    -----
    This assumes the stored interval set and the breakpoint coordinates
    are represented in a compatible coordinate system for your use case.
    """
    hits = []

    for start, end, name in intervals.get(chrom, []):
        if start <= pos <= end:
            hits.append(name)

    return sorted(set(hits))


def annotate_junctions(
    junctions,
    genes_bed: str | None = None,
    exons_bed: str | None = None,
    promoters_bed: str | None = None,
    enhancers_bed: str | None = None,
):
    """
    Annotate harmonized junctions with interval overlaps.

    Added fields in j.annotations
    -----------------------------
    genes_side1, genes_side2
    exons_side1, exons_side2
    promoters_side1, promoters_side2
    enhancers_side1, enhancers_side2
    """
    genes = load_bed4(genes_bed) if genes_bed else None
    exons = load_bed4(exons_bed) if exons_bed else None
    promoters = load_bed4(promoters_bed) if promoters_bed else None
    enhancers = load_bed4(enhancers_bed) if enhancers_bed else None

    for j in junctions:
        ann = j.annotations

        if genes:
            ann["genes_side1"] = point_hits(genes, j.chrom1, j.pos1)
            ann["genes_side2"] = point_hits(genes, j.chrom2, j.pos2)

        if exons:
            ann["exons_side1"] = point_hits(exons, j.chrom1, j.pos1)
            ann["exons_side2"] = point_hits(exons, j.chrom2, j.pos2)

        if promoters:
            ann["promoters_side1"] = point_hits(promoters, j.chrom1, j.pos1)
            ann["promoters_side2"] = point_hits(promoters, j.chrom2, j.pos2)

        if enhancers:
            ann["enhancers_side1"] = point_hits(enhancers, j.chrom1, j.pos1)
            ann["enhancers_side2"] = point_hits(enhancers, j.chrom2, j.pos2)

    return junctions


def annotate_clusters(clusters, junction_by_id):
    """
    Aggregate annotation summaries from junctions into cluster-level annotations.

    Added fields in c.annotations
    -----------------------------
    genes
    exons
    promoters
    enhancers
    """
    for c in clusters:
        genes = set()
        exons = set()
        promoters = set()
        enhancers = set()

        for jid in c.junction_ids:
            j = junction_by_id[jid]
            ann = j.annotations

            for key, collector in (
                ("genes_side1", genes),
                ("genes_side2", genes),
                ("exons_side1", exons),
                ("exons_side2", exons),
                ("promoters_side1", promoters),
                ("promoters_side2", promoters),
                ("enhancers_side1", enhancers),
                ("enhancers_side2", enhancers),
            ):
                for item in ann.get(key, []):
                    collector.add(item)

        c.annotations["genes"] = sorted(genes)
        c.annotations["exons"] = sorted(exons)
        c.annotations["promoters"] = sorted(promoters)
        c.annotations["enhancers"] = sorted(enhancers)

    return clusters
