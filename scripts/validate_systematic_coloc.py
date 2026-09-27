"""Structural and numerical validation of the systematic coloc output."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
PRIORS = {1e-7, 1e-6, 1e-5, 1e-4}


def main():
    coverage = pd.read_csv(ROUND / "tables/systematic_coloc_coverage.csv")
    posterior = pd.read_csv(ROUND / "tables/systematic_coloc_posteriors_all_priors.csv")
    interpretation = pd.read_csv(ROUND / "tables/systematic_coloc_interpretation.csv")
    key = ["gwas_id", "region_id", "outcome"]
    checks = {}
    checks["coverage_rows_1185"] = len(coverage) == 1185
    checks["all_regions_runnable"] = (coverage.analysis_status == "RUN").all()
    checks["all_shared_at_least_50"] = (coverage.n_shared >= 50).all()
    checks["coverage_keys_unique"] = not coverage.duplicated(key).any()
    checks["posterior_rows_4740"] = len(posterior) == 4740
    checks["posterior_keys_unique"] = not posterior.duplicated(key + ["p12"]).any()
    checks["four_priors_per_pair"] = (posterior.groupby(key).size() == 4).all()
    checks["prior_set_exact"] = set(np.round(posterior.p12, 12)) == PRIORS
    checks["no_coloc_errors"] = (posterior.status == "PASS").all()
    pp = posterior[["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]]
    checks["posterior_finite_and_bounded"] = np.isfinite(pp.to_numpy()).all() and (
        (pp.to_numpy() >= 0) & (pp.to_numpy() <= 1)).all()
    max_sum_error = float((pp.sum(axis=1) - 1).abs().max())
    checks["posterior_sum_tolerance"] = max_sum_error < 1e-10
    monotonic_failures = []
    for keys, group in posterior.groupby(key, sort=False):
        ordered = group.sort_values("p12")
        if np.any(np.diff(ordered["PP.H4"].to_numpy()) < -1e-12):
            monotonic_failures.append(dict(zip(key, keys)))
    checks["H4_monotonic_with_p12"] = not monotonic_failures
    checks["interpretation_one_default_row_per_pair"] = len(interpretation) == len(coverage)
    checks = {name: bool(value) for name, value in checks.items()}
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "coverage": {
            "rows": len(coverage), "min_shared": int(coverage.n_shared.min()),
            "median_shared": float(coverage.n_shared.median()),
            "max_shared": int(coverage.n_shared.max()),
            "min_exposure_match_fraction": float(coverage.exposure_match_fraction.min()),
            "median_exposure_match_fraction": float(coverage.exposure_match_fraction.median()),
            "ordered_swapped_matches_total": int(coverage.n_ordered_swapped.sum()),
        },
        "posterior": {
            "rows": len(posterior), "max_sum_error": max_sum_error,
            "runs_with_warnings": int(posterior.warnings.fillna("").astype(bool).sum()),
            "warning_interpretation": "coloc check_dataset minimum-P warnings are retained; they identify weak regional signal and are not computational failures.",
            "monotonic_failures": monotonic_failures,
        },
        "scientific_limit": "ABF is a single-causal-variant regional model. Shared-variant support is not proof of MR direction, mediation, mechanism, or independence from sample overlap.",
    }
    (ROUND / "reports/systematic_coloc_validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report), flush=True)
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
