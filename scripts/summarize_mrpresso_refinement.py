"""Summarize 10,000-draw MR-PRESSO refinement and form a final pair audit."""
from pathlib import Path
import json
import math

import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
ROOT = ROUND / "results/mrpresso"
SEEDS = {271828, 314159}


def parse_p(value):
    text = str(value).strip()
    bound = text.startswith("<")
    return float(text[1:] if bound else text), bound


def read_runs(n_sim):
    rows = []
    for folder in sorted(ROOT.glob(f"*__n{n_sim}")):
        g = pd.read_csv(folder / "global_test.csv").iloc[0]
        main = pd.read_csv(folder / "main_results.csv")
        raw = main.loc[main["MR Analysis"].eq("Raw")].iloc[0]
        corrected = main.loc[main["MR Analysis"].eq("Outlier-corrected")]
        gp, bounded = parse_p(g.Pvalue)
        row = dict(exposure=g.exposure, gwas_id=g.gwas_id, outcome=g.outcome,
                   seed=int(g.seed), n_sim=int(g.n_sim), nsnp_input=int(g.nsnp_input),
                   global_RSS=float(g.RSSobs), global_p_reported=str(g.Pvalue),
                   global_p_numeric_bound=gp, global_p_is_upper_bound=bounded,
                   global_significant=gp < .05, raw_beta=float(raw["Causal Estimate"]),
                   raw_se=float(raw.Sd), raw_p=float(raw["P-value"]),
                   corrected_beta=math.nan, corrected_se=math.nan, corrected_p=math.nan,
                   outlier_markers="", distortion_coefficient=math.nan,
                   distortion_p_reported="")
        if len(corrected):
            c = corrected.iloc[0]
            row.update(corrected_beta=float(c["Causal Estimate"]), corrected_se=float(c.Sd),
                       corrected_p=float(c["P-value"]))
        distortion = folder / "distortion_test.csv"
        if distortion.exists():
            d = pd.read_csv(distortion).iloc[0]
            row.update(outlier_markers=str(d.outlier_markers),
                       distortion_coefficient=float(d.distortion_coefficient),
                       distortion_p_reported=str(d.distortion_p))
        rows.append(row)
    return pd.DataFrame(rows)


def pair_summary(frame, source_precision):
    pairs = []
    for keys, d in frame.groupby(["exposure", "gwas_id", "outcome"], sort=True):
        if len(d) != 2 or set(d.seed) != SEEDS:
            raise RuntimeError(f"Incomplete seeds for {keys}: {d.seed.tolist()}")
        def marker_set(value):
            if pd.isna(value) or not str(value).strip():
                return set()
            return set(filter(None, str(value).split(";")))
        sets = [marker_set(x) for x in d.outlier_markers]
        pairs.append(dict(exposure=keys[0], gwas_id=keys[1], outcome=keys[2],
            final_n_sim_each=int(d.n_sim.iloc[0]), source_precision=source_precision,
            global_significant_both=bool(d.global_significant.all()),
            global_classification_consistent=d.global_significant.nunique() == 1,
            outlier_sets_identical=sets[0] == sets[1],
            robust_outlier_intersection=";".join(sorted(sets[0] & sets[1])),
            stochastic_outlier_union=";".join(sorted(sets[0] | sets[1])),
            n_robust_outlier_intersection=len(sets[0] & sets[1]),
            n_stochastic_outlier_union=len(sets[0] | sets[1]),
            corrected_beta_min=float(d.corrected_beta.min()) if d.corrected_beta.notna().any() else math.nan,
            corrected_beta_max=float(d.corrected_beta.max()) if d.corrected_beta.notna().any() else math.nan,
            corrected_p_max=float(d.corrected_p.max()) if d.corrected_p.notna().any() else math.nan,
            final_interpretation=("GLOBAL_PLEIOTROPY_WITH_SEED_STABLE_OUTLIERS" if sets[0] == sets[1]
                else "GLOBAL_PLEIOTROPY_OUTLIER_SET_REMAINS_SEED_SENSITIVE")))
    return pd.DataFrame(pairs)


def main():
    refined_runs = read_runs(10000)
    if len(refined_runs) != 16:
        raise RuntimeError(f"Expected 16 refined jobs, found {len(refined_runs)}")
    refined_pairs = pair_summary(refined_runs, "REFINED_10000")
    refined_runs.to_csv(ROUND / "tables/mrpresso_refined_seed_results.csv", index=False)
    refined_pairs.to_csv(ROUND / "tables/mrpresso_refined_seed_stability.csv", index=False)

    initial_runs = read_runs(3000)
    initial_pairs = pair_summary(initial_runs, "INITIAL_3000")
    selected = set(zip(refined_pairs.exposure, refined_pairs.outcome))
    stable_initial = initial_pairs.loc[
        ~initial_pairs[["exposure", "outcome"]].apply(tuple, axis=1).isin(selected)].copy()
    final = pd.concat([stable_initial, refined_pairs], ignore_index=True).sort_values(["gwas_id", "outcome"])
    if len(final) != 12 or final.duplicated(["exposure", "outcome"]).any():
        raise RuntimeError("Final MR-PRESSO pair table is not the primary 12 family")
    final.to_csv(ROUND / "tables/mrpresso_final_pair_summary.csv", index=False)
    report = dict(status="PASS", refined_jobs=len(refined_runs), refined_pairs=len(refined_pairs),
        final_pairs=len(final), all_global_classification_consistent=bool(final.global_classification_consistent.all()),
        final_seed_stable_outlier_pairs=int(final.outlier_sets_identical.sum()),
        final_seed_sensitive_outlier_pairs=int((~final.outlier_sets_identical).sum()),
        seed_sensitive_pairs=final.loc[~final.outlier_sets_identical, ["exposure", "outcome"]].to_dict("records"),
        interpretation="The global test is primary for horizontal pleiotropy. If exact outlier sets remain seed-sensitive, no single corrected estimate is privileged; the intersection, union, and corrected range are reported.")
    (ROUND / "reports/mrpresso_refinement_validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
