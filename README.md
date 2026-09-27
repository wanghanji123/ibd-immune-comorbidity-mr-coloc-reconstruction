# IBD immune-comorbidity MR and regional colocalization reconstruction

This repository contains the public-safe code and aggregate outputs for a source-audited reconstruction of genetic-liability associations between inflammatory bowel disease (IBD), Crohn's disease (CD), or ulcerative colitis (UC) and FinnGen R12 ankylosing spondylitis, iridocyclitis, psoriasis, and psoriatic-arthropathy endpoints.

## Interpretive boundary

This is a post hoc methodological reconstruction, not preregistration or independent replication. The results describe directional genetic-liability association and regional shared-variant evidence. They do **not** establish a direct causal effect of diagnosed IBD, mediation by a named gene, or risk of an extraintestinal complication among patients already diagnosed with IBD.

All 12 primary comparisons retained heterogeneity and significant MR-PRESSO global tests. The final evidence gate permits causal confirmation for 0 pairs. High-precision MR-PRESSO retained 6 seed-sensitive exact outlier sets; no preferred seed was selected.

## Data sources

- Exposures: NHGRI-EBI GWAS Catalog GCST004131, GCST004132, and GCST004133 (GRCh37).
- Outcomes: frozen FinnGen R12 endpoints M13_ANKYLOSPON, H7_IRIDOCYCLITIS, L12_PSORIASIS, M13_PSORIARTH, and M13_ANKYLOSPON_STRICT (GRCh38).
- LD reference: audited 1000 Genomes European panel, 502 individuals, MAF > 0.01. This is not FinnGen in-sample LD.

Provider GWAS files and SNP-level effect-size subsets are not redistributed. Obtain them from their original providers after satisfying their access and acknowledgement requirements. The scripts regenerate analysis inputs from those files.

## Repository contents

- `scripts/`: source auditing, instrument construction, liftover, harmonization, MR, MR-PRESSO, systematic colocalization, validation, evidence integration, figures, and supplement generation.
- `config/`: immutable round03 analytical lock.
- `metadata/`: source and software provenance without registration screenshots, email images, credentials, or raw provider data.
- `tables/`: aggregate analysis and validation outputs only.
- `reports/`: machine-readable validation reports and interpretive audits.
- `environment/`: R, Python, PLINK, and Pandoc version records.
- `figures/` and `supplement/`: submission-facing generated artifacts.

## Reproduction sequence

The pipeline uses local paths because provider files are not redistributed. Review and adapt those paths before execution.

1. Audit source files and prepare coordinate chains.
2. Build exposure instruments and independently validate clumping.
3. Reciprocally lift instruments and regions from GRCh37 to GRCh38.
4. Extract exact FinnGen outcome rows and harmonize alleles without proxies.
5. Run frozen MR, robustness analyses, and two-seed MR-PRESSO.
6. Extract dense regional data and run systematic coloc.abf prior sensitivity.
7. Build the conservative evidence gate, figures, and supplement.

See `config/analysis_lock.json` for the exact design and `reports/` for validation outcomes. Do not edit the lock and present the result as preregistered.

## Key limitations

- Exact exposure-outcome participant overlap is unquantified.
- FinnGen outcomes share participants and controls.
- Exposure EAF is unavailable; palindromic variants were excluded and strict binary Steiger analysis was not performed.
- Exposure publication-stage INFO/MAF filtering cannot be completely reconstructed from deposited fields.
- Coloc.abf uses a single-causal-variant model, and the regional threshold is descriptive rather than multiplicity-calibrated across all tests.
- Strict AS is a same-cohort phenotype sensitivity, not independent replication.

## License

Code authored for this reconstruction is released under the MIT License. This license does not apply to third-party GWAS data, software, publications, or provider documentation. Their original terms remain controlling.
