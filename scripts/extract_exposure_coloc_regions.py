"""Stream all prespecified exposure regions from the immutable GRCh37 sources.

This step does not read outcome data and does not select regions by association
strength.  It materialises dense beta/SE records for every region assigned to
each exposure by the frozen primary instrument plan.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from bisect import bisect_right
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
REBUILD = ROUND.parent
RAW = REBUILD / "raw"
RUNTIME = Path("C:/Users/Lenovo/IBD_rebuild_runtime/round03_coloc")
OUT = RUNTIME / "exposure_regions_GRCh37.csv.gz"

SOURCES = {
    "GCST004131": ("IBD", RAW / "ibd_build37_59957_20161107.txt.gz"),
    "GCST004132": ("CD", RAW / "cd_build37_40266_20161107.txt.gz"),
    "GCST004133": ("UC", RAW / "uc_build37_45975_20161107.txt.gz"),
}


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def assigned_regions(plan: pd.DataFrame, gwas_id: str):
    rows = []
    needle = gwas_id + ":"
    for row in plan.itertuples(index=False):
        if any(x.startswith(needle) for x in str(row.index_sources).split(";")):
            rows.append((int(row.chrom37), int(row.start37), int(row.end37), row.region_id))
    return rows


def parse_marker(marker: str):
    try:
        left, ref, alt = marker.split("_")
        chrom, pos = left.split(":")
        chrom = int(chrom)
        pos = int(pos)
        ref = ref.upper()
        alt = alt.upper()
    except (ValueError, AttributeError):
        return None
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        return None
    return chrom, pos, ref, alt


def main():
    RUNTIME.mkdir(parents=True, exist_ok=True)
    plan = pd.read_csv(ROUND / "tables/systematic_coloc_region_plan_GRCh37.csv")
    headers = [
        "trait", "gwas_id", "region_id", "chrom37", "pos37", "ref37", "alt37",
        "marker", "effect_allele37", "other_allele37", "beta", "se", "p",
    ]
    audit = []
    seen = set()
    with gzip.open(OUT, "wt", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=headers)
        writer.writeheader()
        for gwas_id, (trait, source) in SOURCES.items():
            regions = assigned_regions(plan, gwas_id)
            by_chrom = {}
            for chrom, start, end, region_id in regions:
                by_chrom.setdefault(chrom, []).append((start, end, region_id))
            for chrom in by_chrom:
                by_chrom[chrom].sort()
            starts = {chrom: [x[0] for x in values] for chrom, values in by_chrom.items()}
            counts = dict(total_rows=0, malformed_marker=0, outside_assigned_regions=0,
                          invalid_effect=0, duplicate_exact_record=0, retained=0)
            region_counts = {region_id: 0 for *_, region_id in regions}
            with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {"MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P.value"}
                if not required.issubset(reader.fieldnames or []):
                    raise RuntimeError(f"Missing columns in {source}: {required - set(reader.fieldnames or [])}")
                for row in reader:
                    counts["total_rows"] += 1
                    parsed = parse_marker(row["MarkerName"])
                    if parsed is None:
                        counts["malformed_marker"] += 1
                        continue
                    chrom, pos, ref, alt = parsed
                    candidates = by_chrom.get(chrom)
                    if not candidates:
                        counts["outside_assigned_regions"] += 1
                        continue
                    idx = bisect_right(starts[chrom], pos) - 1
                    if idx < 0 or pos > candidates[idx][1]:
                        counts["outside_assigned_regions"] += 1
                        continue
                    region_id = candidates[idx][2]
                    try:
                        beta = float(row["Effect"])
                        se = float(row["StdErr"])
                        p = float(row["P.value"])
                    except (TypeError, ValueError):
                        counts["invalid_effect"] += 1
                        continue
                    if not (math.isfinite(beta) and math.isfinite(se) and se > 0 and
                            math.isfinite(p) and 0 <= p <= 1):
                        counts["invalid_effect"] += 1
                        continue
                    a1 = str(row["Allele1"]).upper()
                    a2 = str(row["Allele2"]).upper()
                    if {a1, a2} != {ref, alt}:
                        counts["invalid_effect"] += 1
                        continue
                    key = (gwas_id, region_id, chrom, pos, ref, alt)
                    if key in seen:
                        counts["duplicate_exact_record"] += 1
                        continue
                    seen.add(key)
                    writer.writerow(dict(trait=trait, gwas_id=gwas_id, region_id=region_id,
                        chrom37=chrom, pos37=pos, ref37=ref, alt37=alt, marker=row["MarkerName"],
                        effect_allele37=a2, other_allele37=a1, beta=beta, se=se, p=p))
                    counts["retained"] += 1
                    region_counts[region_id] += 1
            audit.append(dict(gwas_id=gwas_id, trait=trait, source=str(source),
                source_sha256=sha256(source), assigned_regions=len(regions), counts=counts,
                region_min=min(region_counts.values()), region_max=max(region_counts.values()),
                zero_variant_regions=[k for k, v in region_counts.items() if v == 0]))
            print(json.dumps(audit[-1]), flush=True)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), output=str(OUT),
        output_sha256=sha256(OUT), selection="All finite biallelic SNVs in every prespecified exposure-assigned GRCh37 region; no outcome data read",
        exposures=audit)
    (ROUND / "reports/exposure_coloc_region_extraction.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT), "sha256": report["output_sha256"]}), flush=True)


if __name__ == "__main__":
    main()
