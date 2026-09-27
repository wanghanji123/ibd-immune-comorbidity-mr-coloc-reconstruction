"""Add region, top-variant, and MR context without changing coloc results."""
from pathlib import Path
import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
RUNTIME = Path("C:/Users/Lenovo/IBD_rebuild_runtime/round03_coloc")


def canonical(chrom, pos, ref, alt):
    a, b = sorted((str(ref), str(alt)))
    return f"{chrom}:{int(pos)}:{a}/{b}"


def main():
    result = pd.read_csv(ROUND / "tables/systematic_coloc_interpretation.csv")
    coverage = pd.read_csv(ROUND / "tables/systematic_coloc_coverage.csv")
    plan = pd.read_csv(ROUND / "tables/systematic_coloc_region_plan_GRCh37.csv")
    region38 = pd.read_csv(ROUND / "tables/coloc_region_mapping_GRCh38.csv")
    mr = pd.read_csv(ROUND / "tables/mr_robustness_summary_prePRESSO.csv")
    annotations = []
    for outcome, group in result.groupby("outcome"):
        targets = set(group.top_shared_variant.dropna())
        path = RUNTIME / f"outcome_coloc_{outcome}.csv.gz"
        usecols = ["chrom38", "pos38", "ref38", "alt38", "rsids", "nearest_genes", "extraction_status"]
        for chunk in pd.read_csv(path, usecols=usecols, chunksize=200_000):
            chunk = chunk.loc[chunk.extraction_status.eq("PASS_EXACT_OUTCOME_VARIANT")].copy()
            chunk["top_shared_variant"] = [canonical(*x) for x in
                chunk[["chrom38", "pos38", "ref38", "alt38"]].itertuples(index=False, name=None)]
            hit = chunk.loc[chunk.top_shared_variant.isin(targets),
                ["top_shared_variant", "rsids", "nearest_genes"]].drop_duplicates("top_shared_variant")
            hit["outcome"] = outcome
            annotations.append(hit)
    annotation = pd.concat(annotations, ignore_index=True)
    result = result.merge(annotation, on=["outcome", "top_shared_variant"], how="left", validate="many_to_one")
    result = result.merge(coverage, on=["trait", "gwas_id", "region_id", "outcome", "analysis_role"],
                          how="left", validate="one_to_one", suffixes=("", "_coverage"))
    result = result.merge(plan, on="region_id", how="left", validate="many_to_one")
    result = result.merge(region38, on="region_id", how="left", validate="many_to_one")
    result = result.merge(mr, left_on=["trait", "outcome"], right_on=["exposure", "outcome"],
                          how="left", validate="many_to_one")
    result["regional_signal_context"] = "NEITHER_TRAIT_REACHES_P_LT_1E-5_IN_REGION"
    result.loc[(result.min_p_exposure < 1e-5) | (result.min_p_outcome < 1e-5),
               "regional_signal_context"] = "ONE_TRAIT_REACHES_P_LT_1E-5_IN_REGION"
    result.loc[(result.min_p_exposure < 1e-5) & (result.min_p_outcome < 1e-5),
               "regional_signal_context"] = "BOTH_TRAITS_REACH_P_LT_1E-5_IN_REGION"
    coloc_core = result.evidence_label.eq("PRIOR_ROBUST_STRONG_SHARED_VARIANT_SUPPORT")
    result["coloc_claim_grade"] = result.evidence_label
    weak_trait = ~result.regional_signal_context.eq("BOTH_TRAITS_REACH_P_LT_1E-5_IN_REGION")
    result.loc[coloc_core & weak_trait,
               "coloc_claim_grade"] = "ROBUST_H4_BUT_AT_LEAST_ONE_TRAIT_HAS_WEAK_REGIONAL_SIGNAL"
    mr_core = (result.primary_FDR_positive.fillna(False).astype(bool) &
               result.weighted_median_supports.fillna(False).astype(bool) &
               result.nonMHC_supports.fillna(False).astype(bool) &
               ~result.directional_pleiotropy_signal.fillna(False).astype(bool))
    result["integrated_context_prePRESSO"] = "NO_CONVERGENT_STRONG_CLASSIFICATION"
    result.loc[coloc_core & ~mr_core, "integrated_context_prePRESSO"] = "LOCUS_SHARING_WITHOUT_CORE_MR_SUPPORT"
    result.loc[coloc_core & mr_core, "integrated_context_prePRESSO"] = "CONVERGENT_MR_AND_LOCUS_SUPPORT_WITH_HETEROGENEITY"
    result.to_csv(ROUND / "tables/systematic_coloc_interpretation_annotated.csv", index=False)
    strong = result.loc[coloc_core].copy()
    strong.to_csv(ROUND / "tables/systematic_coloc_prior_robust_strong_regions.csv", index=False)
    print({"annotated_rows": len(result), "top_annotations_missing": int(result.rsids.isna().sum()),
           "robust_strong_rows": len(strong),
           "convergent_prePRESSO": int(strong.integrated_context_prePRESSO.str.startswith("CONVERGENT").sum())})


if __name__ == "__main__":
    main()
