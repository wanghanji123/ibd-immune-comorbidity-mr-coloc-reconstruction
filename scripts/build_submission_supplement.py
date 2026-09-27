"""Build a compact, reviewable Excel supplement from final round03 tables."""
from pathlib import Path
import json

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


ROUND = Path(__file__).resolve().parents[1]
TABLES = ROUND / "tables"
OUT = ROUND / "manuscript_rebuild_2026-09-28/supplement"
OUT.mkdir(parents=True, exist_ok=True)

SHEETS = [
    ("S1_sources", "final_gwas_source_overlap_audit.csv"),
    ("S2_instruments", "instrument_rebuild_summary.csv"),
    ("S3_harmonisation", "harmonisation_summary.csv"),
    ("S4_primary_MR", "mr_primary12.csv"),
    ("S5_MR_robustness", "mr_robustness_summary_prePRESSO.csv"),
    ("S6_PRESSO_runs", "mrpresso_refined_seed_results.csv"),
    ("S7_PRESSO_pairs", "mrpresso_final_pair_summary.csv"),
    ("S8_coloc_loci", "systematic_coloc_prior_robust_strong_regions.csv"),
    ("S9_evidence_gate", "final_integrated_evidence_gate.csv"),
]


def write_frame(ws, frame):
    for j, col in enumerate(frame.columns, 1):
        cell = ws.cell(1, j, col)
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="355C7D")
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    for i, row in enumerate(frame.itertuples(index=False, name=None), 2):
        for j, value in enumerate(row, 1):
            if pd.isna(value):
                value = None
            cell = ws.cell(i, j, value)
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    for j, col in enumerate(frame.columns, 1):
        sample = [str(col)] + [str(x) for x in frame.iloc[:200, j-1].dropna()]
        width = min(42, max(10, max(map(len, sample)) + 2))
        ws.column_dimensions[get_column_letter(j)].width = width


def main():
    missing = [name for _, name in SHEETS if not (TABLES / name).exists()]
    if missing:
        raise RuntimeError("Final tables are not complete: " + ", ".join(missing))

    wb = Workbook()
    readme = wb.active
    readme.title = "README"
    notes = [
        ("Supplement scope", "Compact submission-facing summary. Full SNP and regional posterior files are retained separately with hashes."),
        ("Primary family", "12 IVW comparisons: IBD, CD, and UC x AS, iridocyclitis, psoriasis, and psoriatic arthropathy."),
        ("Causal language", "No row is labelled causal confirmation. Heterogeneity and MR-PRESSO global pleiotropy remain material."),
        ("Colocalization", "Shared-locus support does not establish direction or mediation and uses a single-causal-variant ABF model."),
        ("Sample overlap", "Exact exposure-outcome overlap is unquantified. FinnGen outcomes share a cohort and controls."),
        ("Strict AS", "Same-cohort phenotype sensitivity, not independent replication."),
        ("SNP-level files", "Harmonized inputs and traces are retained in the local audit archive but are not redistributed publicly. Original summary statistics must be obtained from the named providers; the public code regenerates derived inputs."),
        ("Regional outputs", "Full coverage and posterior tables are included in the public aggregate-results archive because they do not redistribute provider effect-size rows."),
    ]
    readme.append(["Item", "Statement"])
    for row in notes:
        readme.append(row)
    write_frame(readme, pd.DataFrame(notes, columns=["Item", "Statement"]))

    counts = {}
    for sheet, filename in SHEETS:
        frame = pd.read_csv(TABLES / filename)
        ws = wb.create_sheet(sheet)
        write_frame(ws, frame)
        counts[sheet] = len(frame)

    path = OUT / "Supplementary_Data_round03_reconstruction.xlsx"
    wb.save(path)
    report = {"status": "PASS", "path": str(path), "sheet_rows": counts}
    (ROUND / "reports/submission_supplement_validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
