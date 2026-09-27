"""Generate submission figures from the final locked round03 tables."""
from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
TABLES = ROUND / "tables"
OUT = ROUND / "manuscript_rebuild_2026-09-28/figures"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 8.5,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

OUTCOME_LABEL = {
    "M13_ANKYLOSPON": "Ankylosing spondylitis",
    "H7_IRIDOCYCLITIS": "Iridocyclitis",
    "L12_PSORIASIS": "Psoriasis",
    "M13_PSORIARTH": "Psoriatic arthropathy",
}
OUTCOME_COLOR = {
    "M13_ANKYLOSPON": "#31708E",
    "H7_IRIDOCYCLITIS": "#2A9D8F",
    "L12_PSORIASIS": "#D97A42",
    "M13_PSORIARTH": "#8E5A9B",
}


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure1_forest(gate):
    order = [(e, o) for e in ["IBD", "CD", "UC"] for o in OUTCOME_LABEL]
    lookup = gate.set_index(["exposure", "outcome"])
    rows = [lookup.loc[x] for x in order]
    fig, (ax, tx) = plt.subplots(1, 2, figsize=(8.1, 6.0),
                                gridspec_kw={"width_ratios": [1.55, 1.05], "wspace": .04})
    y = np.arange(len(rows))[::-1]
    for yi, row, key in zip(y, rows, order):
        color = OUTCOME_COLOR[key[1]]
        ax.errorbar(row.primary_OR, yi,
                    xerr=[[row.primary_OR-row.primary_CI_low],
                          [row.primary_CI_high-row.primary_OR]],
                    fmt="o", color=color, ecolor=color, markersize=5,
                    elinewidth=1.25, capsize=2.2, zorder=3)
    ax.axvline(1, color="#555555", lw=.9, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(.88, 1.43)
    ax.set_xticks([.9, 1.0, 1.1, 1.2, 1.3, 1.4])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    labels = [f"{e}  |  {OUTCOME_LABEL[o]}" for e, o in order]
    ax.set_yticks(y, labels)
    ax.set_xlabel("Odds ratio per 1-unit increase in genetically predicted exposure log odds")
    ax.grid(axis="x", color="#E6E6E6", lw=.7)
    ax.set_title("a  Primary multiplicative random-effects IVW estimates",
                 loc="left", weight="bold", y=1.02)

    tx.set_xlim(0, 1); tx.set_ylim(-.7, len(rows)-.3); tx.axis("off")
    tx.text(.00, len(rows)-.55, "OR (95% CI)", weight="bold", va="bottom")
    tx.text(.60, len(rows)-.55, "BH q", weight="bold", va="bottom")
    for yi, row in zip(y, rows):
        tx.text(.00, yi, f"{row.primary_OR:.3f} ({row.primary_CI_low:.3f}-{row.primary_CI_high:.3f})", va="center")
        tx.text(.60, yi, f"{row.primary_FDR:.2g}", va="center")
    tx.set_title("b  Effect and multiplicity-adjusted evidence",
                 loc="left", weight="bold", y=1.02)
    fig.suptitle("Associations of IBD-family genetic liabilities with FinnGen R12 immune phenotypes",
                 x=.02, ha="left", y=.995, fontsize=11.5, weight="bold")
    fig.text(.02, .01, "Points and horizontal bars are odds ratios and 95% confidence intervals. BH correction covers the locked 12-test primary family.\nAll comparisons showed Cochran-Q heterogeneity and significant MR-PRESSO global tests; statistical association is not causal proof.",
             fontsize=7.5, color="#444444")
    save(fig, "Figure1_primary_MR_forest")


def figure2_evidence_matrix(gate):
    exposures = ["IBD", "CD", "UC"]
    outcomes = list(OUTCOME_LABEL)
    g = gate.set_index(["exposure", "outcome"])
    fill = {
        "MULTIMETHOD_DIRECTIONAL_ASSOCIATION_WITH_INVALIDITY_WARNINGS": "#BFD7EA",
        "METHOD_SENSITIVE_DIRECTIONAL_ASSOCIATION": "#F4C6A8",
        "MHC_SENSITIVE_PRIMARY_NULL": "#F4DFA8",
        "NO_PRIMARY_FDR_SUPPORT": "#E8E8E8",
    }
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for i, e in enumerate(exposures):
        for j, o in enumerate(outcomes):
            r = g.loc[(e, o)]
            rect = plt.Rectangle((j, len(exposures)-1-i), 1, 1,
                                 facecolor=fill[r.directional_MR_grade],
                                 edgecolor="white", lw=2)
            ax.add_patch(rect)
            nloc = int(r.n_prior_robust_shared_regions_with_adequate_signal)
            ax.text(j+.5, len(exposures)-1-i+.61,
                    "Directional MR" if r.primary_FDR < .05 else "Primary MR null",
                    ha="center", va="center", fontsize=8, weight="bold")
            ax.text(j+.5, len(exposures)-1-i+.35,
                    f"Shared loci: {nloc}", ha="center", va="center", fontsize=7.5)
            if r.outlier_sets_identical:
                ax.plot(j+.88, len(exposures)-1-i+.86, marker="o", ms=4,
                        color="#3C3C3C")
            else:
                ax.plot(j+.88, len(exposures)-1-i+.86, marker="o", ms=4,
                        mfc="white", mec="#3C3C3C")
    ax.set_xlim(0, len(outcomes)); ax.set_ylim(0, len(exposures)); ax.set_aspect("equal")
    matrix_labels = ["Ankylosing\nspondylitis", "Iridocyclitis",
                     "Psoriasis", "Psoriatic\narthropathy"]
    ax.set_xticks(np.arange(4)+.5, matrix_labels, rotation=0, ha="center")
    ax.set_yticks(np.arange(3)+.5, exposures[::-1])
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
    ax.set_title("Integrated evidence matrix", loc="left", fontsize=11, weight="bold")
    legend = [
        mpl.patches.Patch(color=fill["MULTIMETHOD_DIRECTIONAL_ASSOCIATION_WITH_INVALIDITY_WARNINGS"], label="Multimethod directional association"),
        mpl.patches.Patch(color=fill["METHOD_SENSITIVE_DIRECTIONAL_ASSOCIATION"], label="Method-sensitive association"),
        mpl.patches.Patch(color=fill["MHC_SENSITIVE_PRIMARY_NULL"], label="Primary null; MHC-sensitive"),
        mpl.patches.Patch(color=fill["NO_PRIMARY_FDR_SUPPORT"], label="No primary FDR support"),
        Line2D([], [], marker="o", color="#3C3C3C", lw=0, label="MR-PRESSO outliers seed-stable"),
        Line2D([], [], marker="o", mfc="white", mec="#3C3C3C", lw=0, label="MR-PRESSO outliers seed-sensitive"),
    ]
    ax.text(0, -.20, "Shared loci require PP.H4 >= 0.80 under both p12=10^-5 and p12=10^-6 and regional P < 10^-5 for both traits.\nEvery cell retains pleiotropy/heterogeneity limitations; no cell denotes causal confirmation.",
            transform=ax.transAxes, fontsize=7.2, color="#444444", va="top")
    ax.legend(handles=legend, ncol=2, frameon=False, loc="upper left",
              bbox_to_anchor=(0, -.35), fontsize=7.2)
    save(fig, "Figure2_integrated_evidence_matrix")


def figure3_coloc():
    d = pd.read_csv(TABLES / "systematic_coloc_prior_robust_strong_regions.csv")
    d = d.loc[(d.analysis_role == "PRIMARY") &
              d.regional_signal_context.eq("BOTH_TRAITS_REACH_P_LT_1E-5_IN_REGION")].copy()
    d["label"] = d.trait + " | " + d.outcome.map(OUTCOME_LABEL) + " | " + d.nearest_genes.fillna("")
    d = d.sort_values(["outcome", "trait", "PP.H4"], ascending=[True, True, False]).reset_index(drop=True)
    h = max(5.0, .24*len(d)+1.4)
    fig, ax = plt.subplots(figsize=(8.0, h))
    y = np.arange(len(d))[::-1]
    ax.hlines(y, d["PP.H4_conservative"], d["PP.H4"], color="#B9B9B9", lw=1.1)
    ax.scatter(d["PP.H4"], y, c="#235789", s=28, label="Default p12=10^-5", zorder=3)
    ax.scatter(d["PP.H4_conservative"], y, c="#D95D39", s=26, marker="D",
               label="Conservative p12=10^-6", zorder=3)
    ax.axvline(.8, color="#555555", ls="--", lw=.8)
    ax.set_xlim(.75, 1.005); ax.set_xlabel("Posterior probability of a shared variant (PP.H4)")
    ax.set_yticks(y, d.label)
    ax.grid(axis="x", color="#E5E5E5", lw=.7)
    ax.legend(frameon=False, ncol=2, loc="lower right")
    ax.set_title("Prior-robust regional shared-variant support", loc="left", fontsize=11, weight="bold")
    fig.text(.02, .006, "Only regions retaining PP.H4 >= 0.80 under both displayed priors and with regional P < 10^-5 in both traits are shown.\nGene labels are nearest-gene annotations for the top shared variant and do not identify a causal gene. ABF assumes at most one causal variant per trait per region.",
             fontsize=7.4, color="#444444")
    save(fig, "Figure3_prior_robust_coloc_loci")


def figure4_workflow():
    fig, ax = plt.subplots(figsize=(8.1, 4.4)); ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
    boxes = [
        (0.35, 3.9, 2.05, 1.15, "Exposure GWAS", "IBD / CD / UC\nGCST004131-3\nGRCh37"),
        (0.35, 1.55, 2.05, 1.35, "FinnGen R12 outcomes", "AS / iridocyclitis /\npsoriasis / PsA\nGRCh38"),
        (3.05, 3.75, 2.25, 1.45, "Instrument pipeline", "Source QC -> GWS ->\nPLINK clumping ->\nreciprocal liftover"),
        (3.05, 1.45, 2.25, 1.55, "Regional pipeline", "130 pre-defined regions ->\ndense extraction -> exact\nallele-pair matching"),
        (5.95, 3.75, 2.15, 1.45, "Directional analyses", "IVW + weighted median +\nEgger + MHC exclusion +\nLOO + MR-PRESSO"),
        (5.95, 1.45, 2.15, 1.55, "Shared-locus analyses", "1,185 region pairs ->\nABF across four priors ->\ncoverage/signal gate"),
        (8.45, 2.55, 1.35, 1.5, "Evidence gate", "Association\nShared locus\nLimitations"),
    ]
    colors = ["#DCEAF4", "#DDEFE9", "#E8EDF3", "#E9F1EC", "#F4E4D9", "#EFE5F1", "#F2F2F2"]
    for (x,y,w,h,title,body), c in zip(boxes, colors):
        ax.add_patch(FancyBboxPatch((x,y), w,h, boxstyle="round,pad=.03,rounding_size=.08",
                                    fc=c, ec="#7A7A7A", lw=.8))
        ax.text(x+w/2, y+h-.25, title, ha="center", va="top", weight="bold", fontsize=8.5)
        ax.text(x+w/2, y+h-.58, body, ha="center", va="top", fontsize=7.4, linespacing=1.25)
    arrows = [((2.4,4.48),(3.05,4.48)),((2.4,2.22),(3.05,2.22)),
              ((5.3,4.48),(5.95,4.48)),((5.3,2.22),(5.95,2.22)),
              ((8.1,4.48),(8.43,3.55)),((8.1,2.22),(8.43,3.05))]
    for a,b in arrows:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", color="#555555", lw=1.0))
    ax.text(.35, 5.65, "Locked reconstruction and evidence-integration workflow", fontsize=12, weight="bold")
    ax.text(.35, .45, "Guardrails: no proxy SNPs; palindromic exposure variants excluded because exposure EAF was unavailable; exact exposure-outcome overlap unquantified;\nstrict AS is same-cohort phenotype sensitivity; coloc is locus-level support and does not validate the MR exclusion restriction.",
            fontsize=7.6, color="#444444")
    save(fig, "Figure4_reconstruction_workflow")


def main():
    gate = pd.read_csv(TABLES / "final_integrated_evidence_gate.csv")
    figure1_forest(gate)
    figure2_evidence_matrix(gate)
    figure3_coloc()
    figure4_workflow()
    print(f"Wrote figures to {OUT}")


if __name__ == "__main__":
    main()
