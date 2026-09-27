options(warn=1)
stopifnot(requireNamespace("data.table", quietly=TRUE))
stopifnot(requireNamespace("coloc", quietly=TRUE))

runtime <- "C:/Users/Lenovo/IBD_rebuild_runtime/round03_coloc"
exposure <- data.table::fread(file.path(runtime, "exposure_coloc_variants_GRCh38.csv.gz"))
exposure[, allele_pair := paste(pmin(ref38, alt38), pmax(ref38, alt38), sep="/")]
exposure[, match_key := paste(chrom38, pos38, allele_pair, sep=":")]
exposure[, duplicate_match_key := duplicated(paste(gwas_id, region_id, match_key, sep="|")) |
  duplicated(paste(gwas_id, region_id, match_key, sep="|"), fromLast=TRUE)]
exposure <- exposure[duplicate_match_key == FALSE]

outcomes <- data.table::data.table(
  outcome=c("M13_ANKYLOSPON", "H7_IRIDOCYCLITIS", "L12_PSORIASIS", "M13_PSORIARTH",
            "M13_ANKYLOSPON_STRICT"),
  analysis_role=c(rep("PRIMARY", 4), "PHENOTYPE_SENSITIVITY_SAME_COHORT"))
priors <- c(1e-5, 1e-6, 1e-7, 1e-4)
minimum_shared <- 50L
posterior_rows <- list(); coverage_rows <- list(); top_rows <- list()
ip <- 0L; ic <- 0L; it <- 0L

for (oi in seq_len(nrow(outcomes))) {
  endpoint <- outcomes$outcome[oi]
  role <- outcomes$analysis_role[oi]
  outcome <- data.table::fread(file.path(runtime, paste0("outcome_coloc_", endpoint, ".csv.gz")))
  outcome <- outcome[extraction_status == "PASS_EXACT_OUTCOME_VARIANT"]
  outcome[, allele_pair := paste(pmin(ref38, alt38), pmax(ref38, alt38), sep="/")]
  outcome[, match_key := paste(chrom38, pos38, allele_pair, sep=":")]
  outcome[, duplicate_match_key := duplicated(paste(region_id, match_key, sep="|")) |
    duplicated(paste(region_id, match_key, sep="|"), fromLast=TRUE)]
  outcome <- outcome[duplicate_match_key == FALSE]
  for (gwas in unique(exposure$gwas_id)) {
    exg <- exposure[gwas_id == gwas]
    trait <- data.table::first(exg$trait)
    for (region in unique(exg$region_id)) {
      ex <- exg[region_id == region]
      oy <- outcome[region_id == region]
      shared <- merge(ex, oy, by="match_key", suffixes=c("_exposure", "_outcome"), sort=FALSE)
      # A canonical allele-pair match may be ordered identically or swapped. No
      # sign alignment is needed for coloc ABF, but the orientation is audited.
      shared[, ordered_alleles_identical := ref38_exposure == ref38_outcome & alt38_exposure == alt38_outcome]
      n_shared <- nrow(shared)
      ic <- ic + 1L
      coverage_rows[[ic]] <- data.table::data.table(trait=trait, gwas_id=gwas, region_id=region,
        outcome=endpoint, analysis_role=role, n_exposure=nrow(ex), n_outcome=nrow(oy),
        n_shared=n_shared, exposure_match_fraction=n_shared/nrow(ex),
        outcome_match_fraction=n_shared/nrow(oy),
        n_ordered_identical=sum(shared$ordered_alleles_identical),
        n_ordered_swapped=sum(!shared$ordered_alleles_identical),
        min_p_exposure=if(nrow(ex)) min(ex$p, na.rm=TRUE) else NA_real_,
        min_p_outcome=if(nrow(oy)) min(oy$pval, na.rm=TRUE) else NA_real_,
        analysis_status=if(n_shared >= minimum_shared) "RUN" else "HOLD_SHARED_VARIANTS_LT_50")
      if (n_shared < minimum_shared) next
      d1 <- list(beta=shared$beta_exposure, varbeta=shared$se_exposure^2,
                 snp=shared$match_key, position=shared$pos38_exposure, type="cc")
      d2 <- list(beta=shared$beta_outcome, varbeta=shared$se_outcome^2,
                 snp=shared$match_key, position=shared$pos38_outcome, type="cc")
      for (p12 in priors) {
        warning_text <- character()
        fit <- tryCatch(withCallingHandlers(
          coloc::coloc.abf(d1, d2, p1=1e-4, p2=1e-4, p12=p12),
          warning=function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }),
          error=function(e) e)
        ip <- ip + 1L
        if (inherits(fit, "error")) {
          posterior_rows[[ip]] <- data.table::data.table(trait=trait, gwas_id=gwas,
            region_id=region, outcome=endpoint, analysis_role=role, p1=1e-4, p2=1e-4,
            p12=p12, nsnps=n_shared, status="ERROR", error=conditionMessage(fit), warnings="")
          next
        }
        s <- fit$summary
        results <- data.table::as.data.table(fit$results)
        top <- results[which.max(SNP.PP.H4)]
        ordered <- results[order(-SNP.PP.H4)]
        credible_n <- which(cumsum(ordered$SNP.PP.H4) >= 0.95)[1]
        posterior_rows[[ip]] <- data.table::data.table(trait=trait, gwas_id=gwas,
          region_id=region, outcome=endpoint, analysis_role=role, p1=1e-4, p2=1e-4,
          p12=p12, nsnps=as.integer(s[["nsnps"]]), PP.H0=as.numeric(s[["PP.H0.abf"]]),
          PP.H1=as.numeric(s[["PP.H1.abf"]]), PP.H2=as.numeric(s[["PP.H2.abf"]]),
          PP.H3=as.numeric(s[["PP.H3.abf"]]), PP.H4=as.numeric(s[["PP.H4.abf"]]),
          top_shared_variant=as.character(top$snp), top_SNP_PP_H4=as.numeric(top$SNP.PP.H4),
          credible_set_95_n=as.integer(credible_n), status="PASS", error="",
          warnings=paste(unique(warning_text), collapse=" | "))
        it <- it + 1L
        top_rows[[it]] <- data.table::data.table(trait=trait, gwas_id=gwas, region_id=region,
          outcome=endpoint, analysis_role=role, p12=p12, match_key=results$snp,
          SNP_PP_H4=results$SNP.PP.H4)
      }
    }
  }
  cat(endpoint, "completed\n")
}

coverage <- data.table::rbindlist(coverage_rows, fill=TRUE)
posterior <- data.table::rbindlist(posterior_rows, fill=TRUE)
top <- data.table::rbindlist(top_rows, fill=TRUE)
stopifnot(nrow(coverage) == data.table::uniqueN(exposure[, .(gwas_id, region_id)]) * nrow(outcomes))
data.table::fwrite(coverage, "tables/systematic_coloc_coverage.csv")
data.table::fwrite(posterior, "tables/systematic_coloc_posteriors_all_priors.csv")
data.table::fwrite(top, file.path(runtime, "systematic_coloc_variant_posteriors.csv.gz"), compress="gzip")

default <- posterior[status == "PASS" & p12 == 1e-5]
conservative <- posterior[status == "PASS" & p12 == 1e-6,
  .(trait, gwas_id, region_id, outcome, PP.H4_conservative=PP.H4)]
interpretation <- merge(default, conservative,
  by=c("trait", "gwas_id", "region_id", "outcome"), all.x=TRUE)
interpretation[, evidence_label := data.table::fcase(
  PP.H4 >= 0.8 & PP.H4_conservative >= 0.8, "PRIOR_ROBUST_STRONG_SHARED_VARIANT_SUPPORT",
  PP.H4 >= 0.8, "DEFAULT_PRIOR_STRONG_BUT_NOT_CONSERVATIVE_ROBUST",
  PP.H4 >= 0.5, "MODERATE_SHARED_VARIANT_SUPPORT",
  PP.H3 >= 0.5, "MORE_SUPPORT_FOR_DISTINCT_VARIANTS_THAN_SHARED_VARIANT",
  default="INCONCLUSIVE_OR_WEAK_REGIONAL_EVIDENCE")]
data.table::fwrite(interpretation, "tables/systematic_coloc_interpretation.csv")
capture.output(sessionInfo(), file="metadata/systematic_coloc_R_sessionInfo.txt")
cat("PASS: coverage rows", nrow(coverage), "posterior rows", nrow(posterior),
    "errors", sum(posterior$status == "ERROR"), "\n")
