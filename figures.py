#!/usr/bin/env python3
"""
figures.py

Generate placeholder figure artifacts for Breakdance events and fusions.

These files are intentionally simple text summaries so downstream reporting can
reference stable figure paths even when no plotting stack is installed.
"""

from __future__ import annotations

from pathlib import Path


def _write_placeholder(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("\n".join(lines) + "\n")


def attach_event_figures(sample: str, events, figure_dir: str | Path | None):
    if figure_dir is None:
        return events
    figure_dir = Path(figure_dir)
    for event in events:
        path = figure_dir / f"{sample}.{event.event_id}.txt"
        _write_placeholder(
            path,
            [
                f"sample\t{sample}",
                f"record_type\tevent",
                f"event_id\t{event.event_id}",
                f"event_type\t{event.event_type}",
                f"label\t{event.event_label or event.event_type}",
                f"location\t{event.location_summary}",
            ],
        )
        event.figure_paths = [str(path)]
    return events


def attach_fusion_figures(sample: str, fusions, figure_dir: str | Path | None):
    if figure_dir is None:
        return fusions
    figure_dir = Path(figure_dir)
    for fusion in fusions:
        path = figure_dir / f"{sample}.{fusion.fusion_id}.txt"
        _write_placeholder(
            path,
            [
                f"sample\t{sample}",
                f"record_type\tfusion",
                f"fusion_id\t{fusion.fusion_id}",
                f"fusion_label\t{fusion.fusion_label}",
                f"location\t{fusion.location_summary}",
            ],
        )
        fusion.figure_path = str(path)
    return fusions


def attach_cluster_figures(sample: str, clusters, figure_dir: str | Path | None):
    if figure_dir is None:
        return clusters
    figure_dir = Path(figure_dir)
    for cluster in clusters:
        path = figure_dir / f"{sample}.{cluster.cluster_id}.txt"
        _write_placeholder(
            path,
            [
                f"sample\t{sample}",
                f"record_type\tcluster",
                f"cluster_id\t{cluster.cluster_id}",
                f"junction_count\t{cluster.junction_count}",
                f"chromosomes\t{','.join(cluster.chromosomes)}",
            ],
        )
        cluster.figure_path = str(path)
    return clusters
