"""Build the locked pair-level evidence gate after formal MR and MR-PRESSO.

This script deliberately separates directional MR evidence from regional
colocalisation.  Colocalisation can support a shared variant at a locus but
cannot repair heterogeneity or horizontal pleiotropy in a genome-wide MR.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
TABLES = ROUND / "tables"
REPORTS = ROUND / "reports"

OUTCOME_LABELS = {
    "M13_ANKYLOSPON": "ankylosing spondylitis",
    "H7_IRIDOCYCLITIS": "iridocyclitis",
    "L12_PSORIASIS": "psoriasis",
    "M13_PSORIARTH": "psoriatic arthropathy",
}


def as_bool(series):
    return series.astype(str).str.lower().map({"true": True, "false": False})


def main():
    mr = pd.read_csv(TABLES / "mr_robustness_summary_prePRESSO.csv")
    presso = pd.read_csv(TABLES / "mrpresso_final_pair_summary.csv")
    coloc = pd.read_csv(TABLES / "systematic_coloc_interpretation_annotated.csv")

    default = coloc.loc[(coloc.analysis_role == "PRIMARY") &
                        np.isclose(coloc.p12, 1e-5)].copy()
    default["robust_coloc"] = default.coloc_claim_grade.eq(
        "PRIOR_ROBUST_STRONG_SHARED_VARIANT_SUPPORT")
    default["adequate_regional_signal"] = default.regional_signal_context.eq(
        "BOTH_TRAITS_REACH_P_LT_1E-5_IN_REGION")

    summaries = []
    for (trait, outcome), d in default.groupby(["trait", "outcome"], sort=True):
        robust = d.loc[d.robust_coloc]
        adequate = robust.loc[robust.adequate_regional_signal]
        summaries.append({
            "exposure": trait,
            "outcome": outcome,
            "n_regions_tested": len(d),
            "n_prior_robust_shared_regions": len(robust),
            "n_prior_robust_shared_regions_with_adequate_signal": len(adequate),
            "robust_shared_region_ids": ";".join(adequate.region_id.astype(str)),
            "robust_shared_region_genes": ";".join(adequate.nearest_genes.fillna("")),
            "max_PP_H4_default": d["PP.H4"].max(),
            "max_PP_H4_conservative": d["PP.H4_conservative"].max(),
        })
    coloc_pair = pd.DataFrame(summaries)

    gate = mr.merge(presso, on=["exposure", "outcome"], how="left",
                    validate="one_to_one").merge(
        coloc_pair, on=["exposure", "outcome"], how="left",
        validate="one_to_one")
    if len(gate) != 12 or gate[["exposure", "outcome"]].duplicated().any():
        raise RuntimeError("Evidence gate must contain the locked primary 12 pairs")

    for col in ["primary_FDR_positive", "heterogeneity_present",
                "weighted_median_supports", "nonMHC_supports",
                "directional_pleiotropy_signal"]:
        gate[col] = as_bool(gate[col])

    gate["outcome_label"] = gate.outcome.map(OUTCOME_LABELS)
    gate["directional_MR_grade"] = "NO_PRIMARY_FDR_SUPPORT"
    broad = (gate.primary_FDR_positive & gate.weighted_median_supports &
             gate.nonMHC_supports & ~gate.directional_pleiotropy_signal)
    gate.loc[broad, "directional_MR_grade"] = (
        "MULTIMETHOD_DIRECTIONAL_ASSOCIATION_WITH_INVALIDITY_WARNINGS")
    partial = gate.primary_FDR_positive & ~broad
    gate.loc[partial, "directional_MR_grade"] = (
        "METHOD_SENSITIVE_DIRECTIONAL_ASSOCIATION")

    # A primary-null result that appears only after MHC removal is instability,
    # not upgraded evidence.
    mhc_sensitive = ~gate.primary_FDR_positive & gate.nonMHC_supports
    gate.loc[mhc_sensitive, "directional_MR_grade"] = (
        "MHC_SENSITIVE_PRIMARY_NULL")

    gate["coloc_grade"] = np.where(
        gate.n_prior_robust_shared_regions_with_adequate_signal.fillna(0).gt(0),
        "PRIOR_ROBUST_SHARED_LOCUS_SUPPORT",
        "NO_PRIOR_ROBUST_SHARED_LOCUS_SUPPORT")

    gate["invalidity_flag"] = (
        gate.heterogeneity_present | gate.global_significant_both |
        gate.directional_pleiotropy_signal)
    gate["causal_claim_permitted"] = False

    def integrated(row):
        mr_positive = row.directional_MR_grade == (
            "MULTIMETHOD_DIRECTIONAL_ASSOCIATION_WITH_INVALIDITY_WARNINGS")
        mr_sensitive = row.directional_MR_grade == "METHOD_SENSITIVE_DIRECTIONAL_ASSOCIATION"
        shared = row.coloc_grade == "PRIOR_ROBUST_SHARED_LOCUS_SUPPORT"
        if mr_positive and shared:
            return "CONVERGENT_DIRECTIONAL_AND_SHARED_LOCUS_EVIDENCE_WITH_PLEIOTROPY_LIMITATION"
        if mr_positive:
            return "DIRECTIONAL_ASSOCIATION_WITHOUT_ROBUST_SHARED_LOCUS_SUPPORT"
        if mr_sensitive and shared:
            return "METHOD_SENSITIVE_DIRECTIONAL_SIGNAL_WITH_SHARED_LOCUS_SUPPORT"
        if mr_sensitive:
            return "METHOD_SENSITIVE_DIRECTIONAL_SIGNAL_ONLY"
        if shared:
            return "SHARED_LOCUS_EVIDENCE_WITHOUT_STABLE_DIRECTIONAL_MR"
        if row.directional_MR_grade == "MHC_SENSITIVE_PRIMARY_NULL":
            return "PRIMARY_NULL_BUT_MHC_SENSITIVE"
        return "NO_CONVERGENT_PRIMARY_EVIDENCE"

    gate["integrated_evidence_grade"] = gate.apply(integrated, axis=1)
    gate["reporting_guardrail"] = (
        "Describe genetic-liability association/shared architecture; do not "
        "claim causal proof. MR-PRESSO global pleiotropy and heterogeneity "
        "remain relevant even if outlier-corrected estimates are non-null."
    )

    ordered = [
        "exposure", "gwas_id", "outcome", "outcome_label", "nsnp",
        "primary_OR", "primary_CI_low", "primary_CI_high", "primary_p",
        "primary_FDR", "heterogeneity_Q", "heterogeneity_Q_df",
        "heterogeneity_Q_p", "weighted_median_OR", "weighted_median_p",
        "egger_OR", "egger_p", "egger_intercept", "egger_intercept_p",
        "nonMHC_nsnp", "nonMHC_OR", "nonMHC_p", "strict_AS_OR",
        "strict_AS_p", "loo_beta_min", "loo_beta_max",
        "global_significant_both", "final_n_sim_each", "source_precision",
        "outlier_sets_identical", "robust_outlier_intersection",
        "stochastic_outlier_union", "corrected_beta_min", "corrected_beta_max",
        "corrected_p_max", "n_regions_tested",
        "n_prior_robust_shared_regions",
        "n_prior_robust_shared_regions_with_adequate_signal",
        "robust_shared_region_ids", "robust_shared_region_genes",
        "max_PP_H4_default", "max_PP_H4_conservative",
        "directional_MR_grade", "coloc_grade", "invalidity_flag",
        "causal_claim_permitted", "integrated_evidence_grade",
        "reporting_guardrail",
    ]
    gate[ordered].sort_values(["exposure", "outcome"]).to_csv(
        TABLES / "final_integrated_evidence_gate.csv", index=False)

    report = {
        "status": "PASS",
        "n_locked_primary_pairs": len(gate),
        "n_primary_fdr_positive": int(gate.primary_FDR_positive.sum()),
        "n_global_pleiotropy_positive": int(gate.global_significant_both.sum()),
        "n_heterogeneity_positive": int(gate.heterogeneity_present.sum()),
        "n_pairs_with_prior_robust_shared_locus": int(
            gate.n_prior_robust_shared_regions_with_adequate_signal.fillna(0).gt(0).sum()),
        "n_pairs_permitting_causal_claim": int(gate.causal_claim_permitted.sum()),
        "grade_counts": gate.integrated_evidence_grade.value_counts().to_dict(),
        "interpretive_boundary": (
            "The evidence gate is conservative by design. Colocalisation is "
            "treated as locus-level shared-variant support and never as proof "
            "that the genome-wide MR exclusion restriction holds."
        ),
    }
    (REPORTS / "final_integrated_evidence_gate_validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
