#!/usr/bin/env python3
from __future__ import annotations
import gzip
import re
from pathlib import Path
from typing import Optional
from schema import RawSVCall

def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if chrom.startswith("chr"):
        return chrom
    if chrom in {"X", "Y", "M", "MT"} or chrom.isdigit():
        return f"chr{chrom}"
    return chrom

def normalize_svtype(svtype: str) -> str:
    if not svtype:
        return "NA"
    svtype = svtype.upper()
    if svtype == "TRA":
        return "BND"
    return svtype

def parse_info_field(info_str: str) -> dict[str, str]:
    info = {}
    if not info_str or info_str == ".":
        return info
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = "True"
    return info

def parse_bnd_alt(alt: str) -> tuple[Optional[str], Optional[int], Optional[str], Optional[str]]:
    m = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
    if not m:
        return None, None, None, None
    chrom2 = normalize_chrom(m.group(1))
    pos2 = int(m.group(2))
    strand1 = None
    strand2 = None
    if alt.startswith("[") or alt.startswith("]"):
        strand1 = "-"
    elif alt.endswith("[") or alt.endswith("]"):
        strand1 = "+"
    if "[" in alt:
        strand2 = "+"
    elif "]" in alt:
        strand2 = "-"
    return chrom2, pos2, strand1, strand2

def extract_support(info: dict[str, str]) -> Optional[int]:
    for key in ("RE", "SUPPORT", "SU", "DV", "SR"):
        if key in info:
            raw = info[key].split(",")[0]
            try:
                return int(float(raw))
            except ValueError:
                continue
    return None

def extract_svlen(info: dict[str, str]) -> Optional[int]:
    if "SVLEN" in info:
        raw = info["SVLEN"].split(",")[0]
        try:
            return int(float(raw))
        except ValueError:
            return None
    return None

def _open_textmaybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")

def parse_vcf_record(line: str, caller: str, source_path: Path, sample_name: str):
    if not line.strip() or line.startswith("#"):
        return None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        return None
    chrom1 = normalize_chrom(fields[0])
    pos1 = int(fields[1])
    record_id = fields[2]
    alt = fields[4]
    qual_str = fields[5]
    filt = fields[6]
    info_str = fields[7]
    info = parse_info_field(info_str)
    svtype = normalize_svtype(info.get("SVTYPE", "NA"))
    chrom2 = chrom1
    pos2 = pos1
    strand1 = None
    strand2 = None
    if svtype == "BND":
        c2, p2, s1, s2 = parse_bnd_alt(alt)
        if c2 is not None and p2 is not None:
            chrom2 = c2
            pos2 = p2
            strand1 = s1
            strand2 = s2
    else:
        if "CHR2" in info:
            chrom2 = normalize_chrom(info["CHR2"])
        if "END" in info:
            try:
                pos2 = int(float(info["END"].split(",")[0]))
            except ValueError:
                pos2 = pos1
    qual = None
    if qual_str not in {".", ""}:
        try:
            qual = float(qual_str)
        except ValueError:
            qual = None
    return RawSVCall(
        caller=caller,
        record_id=record_id,
        svtype=svtype,
        chrom1=chrom1,
        pos1=pos1,
        chrom2=chrom2,
        pos2=pos2,
        strand1=strand1,
        strand2=strand2,
        svlen=extract_svlen(info),
        support=extract_support(info),
        qual=qual,
        filt=filt,
        source_path=str(source_path),
        sample_name=sample_name,
        info=info,
    )

def read_vcf_calls(vcf_path: Path, caller: str, sample_name: str):
    calls = []
    with _open_textmaybe_gzip(vcf_path) as handle:
        for line in handle:
            rec = parse_vcf_record(line, caller, vcf_path, sample_name)
            if rec is not None:
                calls.append(rec)
    return calls

def discover_caller_vcfs(sample_output_dir: Path):
    candidates = {
        "sniffles2": [
            sample_output_dir / "sniffles2" / "sniffles2.vcf",
            sample_output_dir / "sniffles2" / "sniffles2.vcf.gz",
        ],
        "cutesv": [
            sample_output_dir / "cutesv" / "cutesv.vcf",
            sample_output_dir / "cutesv" / "cutesv.vcf.gz",
        ],
        "svim": [
            sample_output_dir / "svim" / "variants.vcf",
            sample_output_dir / "svim" / "variants.vcf.gz",
        ],
    }
    found = {}
    for caller, paths in candidates.items():
        for p in paths:
            if p.exists() and p.is_file():
                found[caller] = p
                break
    return found
