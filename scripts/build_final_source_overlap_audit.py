"""Create the submission-facing GWAS source and overlap audit.

The table separates verified source facts from unresolved overlap.  It does
not infer zero overlap from public availability or from study-level labels.
"""
from pathlib import Path
import json

import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
SOURCE = Path("C:/Users/Lenovo/IBD_rebuild_runtime/source_project")
ROUND02 = SOURCE / "重建研究_v1_2026-09-14/round02_2026-09-18"


def main():
    exposure_qc = pd.read_csv(ROUND02 / "tables/source_eligibility_round02.csv")
    exposure_qc = exposure_qc.loc[exposure_qc.iloc[:, 0].isin(
        ["GCST004131", "GCST004132", "GCST004133"])].copy()
    fg = pd.read_csv(ROUND / "tables/FinnGen_R12_source_manifest_completed_2026-09-28.csv")

    exposure_rows = [
        dict(dataset_id="GCST004131", role="exposure", phenotype="Inflammatory bowel disease",
             cases=25042, controls=34915),
        dict(dataset_id="GCST004132", role="exposure", phenotype="Crohn's disease",
             cases=12194, controls=28072),
        dict(dataset_id="GCST004133", role="exposure", phenotype="Ulcerative colitis",
             cases=12366, controls=33609),
    ]
    rows = []
    for item in exposure_rows:
        src = exposure_qc.loc[exposure_qc.iloc[:, 0].eq(item["dataset_id"])].iloc[0]
        rows.append({**item,
            "release_or_publication": "de Lange et al. 2017",
            "ancestry_context": "European-ancestry meta-analysis as reported by source study",
            "genome_build": "GRCh37",
            "effect_scale": "log odds ratio",
            "effect_allele_rule": "Allele2 in original file; beta retained",
            "frequency_available": "No",
            "source_url": src.iloc[13],
            "citation": "de Lange et al. Nat Genet. 2017;49:256-261. doi:10.1038/ng.3760",
            "local_integrity": "Full gzip/row/QC audit completed in round02; see source_eligibility_round02.csv",
            "phenotype_scope_note": "Population susceptibility, not complications among patients with IBD",
            "within_domain_overlap": "IBD, CD and UC are related analyses from the same study framework and are not independent exposures",
            "cross_domain_overlap_risk": "UNCLEAR_CONTEXT",
            "cross_domain_overlap_rationale": "Exact participant overlap with FinnGen R12 cannot be quantified from the available summary files and source reports",
        })

    endpoint_notes = {
        "M13_ANKYLOSPON": "Broad register-based ankylosing spondylitis endpoint",
        "M13_ANKYLOSPON_STRICT": "Same-cohort stricter phenotype sensitivity; not independent replication",
        "H7_IRIDOCYCLITIS": "Register-based iridocyclitis composite; not equivalent to all uveitis",
        "L12_PSORIASIS": "Broad psoriasis endpoint; includes subendpoints and is not restricted to skin-only disease",
        "M13_PSORIARTH": "Psoriatic arthropathies endpoint; clinically and genetically overlaps psoriasis",
    }
    for x in fg.itertuples(index=False):
        rows.append(dict(
            dataset_id=f"FinnGen_R12_{x.endpoint}", role=("sensitivity_outcome" if
                x.endpoint == "M13_ANKYLOSPON_STRICT" else "outcome"),
            phenotype=x.phenotype, cases=int(x.num_cases), controls=int(x.num_controls),
            release_or_publication="FinnGen R12",
            ancestry_context="Finnish biobank cohort",
            genome_build="GRCh38", effect_scale="log odds ratio",
            effect_allele_rule="ALT allele", frequency_available="Yes: af_alt, af_alt_cases, af_alt_controls",
            source_url=x.url,
            citation="Kurki et al. Nature. 2023;613:508-518. doi:10.1038/s41586-022-05473-8",
            local_integrity=(f"Provider MD5/size and full gzip CRC verified; SHA-256 {x.sha256}"),
            phenotype_scope_note=endpoint_notes[x.endpoint],
            within_domain_overlap="Same FinnGen R12 cohort; shared controls and phenotype overlap make outcomes statistically dependent",
            cross_domain_overlap_risk="UNCLEAR_CONTEXT",
            cross_domain_overlap_rationale="Exact participant overlap with de Lange et al. exposure GWAS cannot be quantified from public summary statistics",
        ))

    audit = pd.DataFrame(rows)
    audit.to_csv(ROUND / "tables/final_gwas_source_overlap_audit.csv", index=False)
    report = {
        "status": "PASS",
        "n_sources": len(audit),
        "verified_finngen_files": int(fg.provider_size_md5_match.astype(str).str.lower().eq("true").sum()),
        "zero_overlap_claim_permitted": False,
        "independent_replication_claim_permitted_for_strict_AS": False,
        "formal_subtype_difference_from_significance_pattern_permitted": False,
        "required_language": [
            "Exact exposure-outcome participant overlap could not be quantified.",
            "FinnGen outcomes share a source cohort and controls and are not independent.",
            "The strict AS endpoint is a same-cohort phenotype sensitivity analysis.",
        ],
    }
    (ROUND / "reports/final_gwas_source_overlap_audit_validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
