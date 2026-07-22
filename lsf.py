#!/usr/bin/env python3
"""
lsf.py

Minimal helper for writing an LSF wrapper script around Breakdance execution.
"""

from __future__ import annotations

from pathlib import Path


def write_lsf_launcher(script_path: str | Path, command: str, job_name: str = "breakdance") -> str:
    script_path = Path(script_path)
    script_path.parent.mkdir(parents=True, exist_ok=True)
    content = "\n".join(
        [
            "#!/bin/bash",
            f"#BSUB -J {job_name}",
            "#BSUB -q normal",
            "#BSUB -n 1",
            "#BSUB -R span[hosts=1]",
            "set -euo pipefail",
            command,
            "",
        ]
    )
    script_path.write_text(content)
    return str(script_path)
