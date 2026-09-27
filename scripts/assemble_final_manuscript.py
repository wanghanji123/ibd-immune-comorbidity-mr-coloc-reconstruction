"""Fill result-derived prose and main tables in the manuscript source."""
from pathlib import Path
import math

import numpy as np
import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
M = ROUND / "manuscript_rebuild_2026-09-28"
T = ROUND / "tables"

OUTCOME = {
    "M13_ANKYLOSPON": "AS",
    "H7_IRIDOCYCLITIS": "Iridocyclitis",
    "L12_PSORIASIS": "Psoriasis",
    "M13_PSORIARTH": "Psoriatic arthropathy",
}


def pformat(x):
    if pd.isna(x):
        return "NA"
    if x < .001:
        exponent = int(math.floor(math.log10(x)))
        mantissa = x / 10**exponent
        return f"{mantissa:.2f} x 10^{exponent}"
    return f"{x:.3f}"


def presso_text(presso):
    stable = int(presso.outlier_sets_identical.sum())
    sensitive = len(presso) - stable
    corrected_sig = int((presso.corrected_p_max < .05).sum())
    refined = int((presso.source_precision == "REFINED_10000").sum())
    names = [f"{r.exposure}-{OUTCOME[r.outcome]}" for r in
             presso.loc[~presso.outlier_sets_identical].itertuples(index=False)]
    return (
        f"The MR-PRESSO global test was significant under both fixed seeds for all "
        f"12 primary comparisons. Eight initially seed-unstable pairs were refined "
        f"to 10,000 simulated distributions per seed; the other four retained their "
        f"3,000-distribution runs because their exact outlier sets were already identical. "
        f"After refinement, exact outlier sets were identical for {stable} pairs and "
        f"remained seed-sensitive for {sensitive}: {', '.join(names) if names else 'none'}. "
        f"The maximum outlier-corrected P value across the two retained seeds was below "
        f"0.05 for {corrected_sig} pairs. Corrected estimates were reported as seed-specific "
        f"ranges rather than selecting the more favorable run. Because every global test "
        f"was significant, outlier correction did not restore an assumption of no horizontal "
        f"pleiotropy and did not upgrade any comparison to causal confirmation."
    )


def evidence_text(gate):
    convergent = gate.integrated_evidence_grade.eq(
        "CONVERGENT_DIRECTIONAL_AND_SHARED_LOCUS_EVIDENCE_WITH_PLEIOTROPY_LIMITATION")
    directional_only = gate.integrated_evidence_grade.eq(
        "DIRECTIONAL_ASSOCIATION_WITHOUT_ROBUST_SHARED_LOCUS_SUPPORT")
    shared_only = gate.integrated_evidence_grade.eq(
        "SHARED_LOCUS_EVIDENCE_WITHOUT_STABLE_DIRECTIONAL_MR")
    method_sensitive = gate.integrated_evidence_grade.str.startswith("METHOD_SENSITIVE")
    return (
        f"The conservative evidence gate identified {int(convergent.sum())} pairs with both "
        f"multimethod directional association and prior-robust shared-locus support, "
        f"{int(directional_only.sum())} with directional association but no robust shared "
        f"locus, {int(shared_only.sum())} with shared-locus support but no stable primary "
        f"directional MR, and {int(method_sensitive.sum())} with a method-sensitive primary "
        f"signal. No pair permitted a causal-confirmation label because all retained "
        f"heterogeneity and MR-PRESSO global-test warnings. Pair-specific classifications "
        f"are provided in Table 3 and Supplementary Table S9."
    )


def make_tables(gate, mr, instruments, sources):
    source_lookup = sources.set_index("dataset_id")
    rows = []
    for trait, gwas in [("IBD", "GCST004131"), ("CD", "GCST004132"), ("UC", "GCST004133")]:
        inst = instruments.loc[(instruments.trait == trait) &
                               (instruments.branch == "primary_all")].iloc[0]
        src = source_lookup.loc[gwas]
        rows.append({"Exposure": trait, "GWAS": gwas,
                     "Cases/controls": f"{int(src.cases):,}/{int(src.controls):,}",
                     "Pre-harmonisation instruments": int(inst.clumped_index),
                     "Minimum F": f"{inst.min_F:.2f}",
                     "Median F": f"{inst.median_F:.2f}"})
    t1 = pd.DataFrame(rows).to_markdown(index=False)

    t2 = mr.copy()
    t2["Exposure"] = t2.exposure
    t2["Outcome"] = t2.outcome.map(OUTCOME)
    t2["SNPs"] = t2.nsnp.astype(int)
    t2["IVW OR (95% CI)"] = t2.apply(
        lambda r: f"{r.OR:.3f} ({r.CI_low:.3f}-{r.CI_high:.3f})", axis=1)
    t2["P"] = t2.p.map(pformat)
    t2["BH q"] = t2.FDR_BH_primary12.map(pformat)
    t2["Q P"] = t2.Q_p.map(pformat)
    t2 = t2[["Exposure", "Outcome", "SNPs", "IVW OR (95% CI)", "P", "BH q", "Q P"]]
    t2 = t2.sort_values(["Exposure", "Outcome"]).to_markdown(index=False)

    t3 = gate.copy()
    t3["Exposure"] = t3.exposure
    t3["Outcome"] = t3.outcome.map(OUTCOME)
    t3["Directional MR"] = t3.directional_MR_grade.str.replace("_", " ").str.title()
    t3["Robust shared loci"] = t3.n_prior_robust_shared_regions_with_adequate_signal.astype(int)
    t3["PRESSO outliers"] = np.where(t3.outlier_sets_identical,
                                      "Seed-stable", "Seed-sensitive")
    t3["Integrated interpretation"] = t3.integrated_evidence_grade.str.replace("_", " ").str.title()
    t3 = t3[["Exposure", "Outcome", "Directional MR", "Robust shared loci",
             "PRESSO outliers", "Integrated interpretation"]]
    t3 = t3.sort_values(["Exposure", "Outcome"]).to_markdown(index=False)

    pair_coloc = gate[["exposure", "outcome", "n_regions_tested",
                       "n_prior_robust_shared_regions_with_adequate_signal",
                       "robust_shared_region_genes"]].copy()
    pair_coloc.columns = ["Exposure", "Outcome", "Regions tested",
                          "Prior-robust regions", "Nearest-gene annotations"]
    pair_coloc["Outcome"] = pair_coloc.Outcome.map(OUTCOME)
    pair_coloc["Prior-robust regions"] = pair_coloc["Prior-robust regions"].astype(int)
    pair_coloc = pair_coloc.sort_values(["Exposure", "Outcome"]).to_markdown(index=False)

    return (
        "**Table 1. Exposure sources and primary instrument construction.**\n\n" + t1 +
        "\n\n**Table 2. Locked primary multiplicative random-effects IVW estimates.**\n\n" + t2 +
        "\n\n**Table 3. Conservative integrated evidence gate.**\n\n" + t3 +
        "\n\n**Table 4. Pair-level systematic colocalization summary.**\n\n" + pair_coloc +
        "\n\nTable notes: AS, ankylosing spondylitis; BH, Benjamini-Hochberg; "
        "CD, Crohn's disease; CI, confidence interval; IBD, inflammatory bowel "
        "disease; IVW, inverse-variance weighted; MHC, major histocompatibility "
        "complex; UC, ulcerative colitis. Shared-locus counts require PP.H4 >= "
        "0.80 under p12 = 10^-5 and 10^-6 and regional P < 10^-5 in both traits."
    )


def main():
    source = (M / "manuscript_source.md").read_text(encoding="utf-8")
    presso = pd.read_csv(T / "mrpresso_final_pair_summary.csv")
    gate = pd.read_csv(T / "final_integrated_evidence_gate.csv")
    mr = pd.read_csv(T / "mr_primary12.csv")
    instruments = pd.read_csv(T / "instrument_rebuild_summary.csv")
    sources = pd.read_csv(T / "final_gwas_source_overlap_audit.csv")
    code_text = (
        "Analysis scripts, aggregate outputs, source manifests, validation reports, "
        "and environment records are available at "
        "https://github.com/wanghanji123/ibd-immune-comorbidity-mr-coloc-reconstruction "
        "under an "
        "MIT license. Third-party GWAS summary statistics and SNP-level effect-size "
        "subsets are excluded from the public repository and must be obtained from "
        "their original providers; the documented scripts regenerate the analysis "
        "inputs after those files are placed in the declared local paths."
    )
    replacements = {
        "{{MRPRESSO_FINAL_RESULTS}}": presso_text(presso),
        "{{INTEGRATED_EVIDENCE_GATE_RESULTS}}": evidence_text(gate),
        "{{MAIN_TABLES}}": make_tables(gate, mr, instruments, sources),
        "{{CODE_ARCHIVE_AND_RELEASE}}": code_text,
    }
    for marker, value in replacements.items():
        if marker not in source:
            raise RuntimeError(f"Missing manuscript marker: {marker}")
        source = source.replace(marker, value)
    unresolved = [x for x in ["{{MRPRESSO", "{{INTEGRATED", "{{MAIN_TABLES}}",
                                  "{{CODE_ARCHIVE"] if x in source]
    if unresolved:
        raise RuntimeError(f"Unresolved quantitative markers: {unresolved}")
    (M / "manuscript_assembled.md").write_text(source, encoding="utf-8")
    print({"status": "PASS", "rows_gate": len(gate), "rows_presso": len(presso)})


if __name__ == "__main__":
    main()
