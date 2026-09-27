options(warn=1)
stopifnot(requireNamespace("data.table", quietly=TRUE))
stopifnot(requireNamespace("Rsamtools", quietly=TRUE))
stopifnot(requireNamespace("GenomicRanges", quietly=TRUE))
stopifnot(requireNamespace("IRanges", quietly=TRUE))

project <- "C:/Users/Lenovo/IBD_rebuild_runtime/source_project"
finngen <- "C:/Users/Lenovo/IBD_rebuild_runtime/finngen_R12"
runtime <- "C:/Users/Lenovo/IBD_rebuild_runtime/round03_coloc"
dir.create(runtime, recursive=TRUE, showWarnings=FALSE)
sources <- data.table::data.table(
  outcome=c("M13_ANKYLOSPON", "H7_IRIDOCYCLITIS", "L12_PSORIASIS", "M13_PSORIARTH",
            "M13_ANKYLOSPON_STRICT"),
  analysis_role=c(rep("PRIMARY", 4), "PHENOTYPE_SENSITIVITY_SAME_COHORT"),
  path=c(
    file.path(project, "data_external/finngen_as_r12/finngen_R12_M13_ANKYLOSPON.gz"),
    file.path(finngen, "finngen_R12_H7_IRIDOCYCLITIS.gz"),
    file.path(finngen, "finngen_R12_L12_PSORIASIS.gz"),
    file.path(finngen, "finngen_R12_M13_PSORIARTH.gz"),
    file.path(project, "data_external/finngen_as_r12_strict/finngen_R12_M13_ANKYLOSPON_STRICT.gz")))

exposure <- data.table::fread(file.path(runtime, "exposure_coloc_variants_GRCh38.csv.gz"))
regions <- exposure[, .(n_chrom=data.table::uniqueN(chrom38), chrom38=data.table::first(chrom38),
                        start38=min(pos38), end38=max(pos38), mapped_exposure_variants=.N),
                    by=region_id]
stopifnot(all(regions$n_chrom == 1L))
data.table::setorder(regions, chrom38, start38, end38)
ranges <- GenomicRanges::GRanges(as.character(regions$chrom38),
  IRanges::IRanges(regions$start38, regions$end38), strand="*")
header <- c("chrom38", "pos38", "ref38", "alt38", "rsids", "nearest_genes", "pval",
            "mlogp", "beta", "se", "af_alt", "af_alt_cases", "af_alt_controls")
summary_rows <- list()

for (i in seq_len(nrow(sources))) {
  endpoint <- sources$outcome[i]
  role <- sources$analysis_role[i]
  path <- sources$path[i]
  stopifnot(file.exists(path), file.exists(paste0(path, ".tbi")))
  raw_lines <- Rsamtools::scanTabix(Rsamtools::TabixFile(path), param=ranges)
  stopifnot(length(raw_lines) == nrow(regions))
  lengths_by_region <- lengths(raw_lines)
  nonempty <- which(lengths_by_region > 0L)
  pieces <- vector("list", length(nonempty))
  for (j in seq_along(nonempty)) {
    k <- nonempty[j]
    text <- paste(raw_lines[[k]], collapse="\n")
    x <- data.table::fread(text=text, sep="\t", header=FALSE, col.names=header,
                          showProgress=FALSE)
    x[, region_id := regions$region_id[k]]
    pieces[[j]] <- x
  }
  dat <- data.table::rbindlist(pieces, use.names=TRUE)
  dat[, `:=`(outcome=endpoint, analysis_role=role)]
  dat[, valid_effect := is.finite(beta) & is.finite(se) & se > 0 & is.finite(pval) &
         pval >= 0 & pval <= 1 & nchar(ref38) == 1L & nchar(alt38) == 1L &
         ref38 %in% c("A","C","G","T") & alt38 %in% c("A","C","G","T") &
         ref38 != alt38]
  dat[, variant_key := paste(chrom38, pos38, ref38, alt38, sep=":")]
  dat[, duplicate_region_key := duplicated(paste(region_id, variant_key, sep="|")) |
         duplicated(paste(region_id, variant_key, sep="|"), fromLast=TRUE)]
  dat[, extraction_status := data.table::fcase(
    !valid_effect, "HOLD_INVALID_EFFECT_OR_NON_SNV",
    duplicate_region_key, "HOLD_DUPLICATE_REGION_KEY",
    default="PASS_EXACT_OUTCOME_VARIANT")]
  out_path <- file.path(runtime, paste0("outcome_coloc_", endpoint, ".csv.gz"))
  data.table::fwrite(dat, out_path, compress="gzip")
  summary_rows[[i]] <- dat[, .(queried_regions=nrow(regions), raw_rows=.N,
    pass_rows=sum(extraction_status == "PASS_EXACT_OUTCOME_VARIANT"),
    invalid_rows=sum(extraction_status == "HOLD_INVALID_EFFECT_OR_NON_SNV"),
    duplicate_rows=sum(extraction_status == "HOLD_DUPLICATE_REGION_KEY"),
    empty_regions=sum(lengths_by_region == 0L)), by=.(outcome, analysis_role)]
  cat(endpoint, "raw", nrow(dat), "pass", sum(dat$extraction_status == "PASS_EXACT_OUTCOME_VARIANT"), "\n")
}
summary <- data.table::rbindlist(summary_rows)
data.table::fwrite(regions, "tables/coloc_region_mapping_GRCh38.csv")
data.table::fwrite(summary, "tables/outcome_coloc_region_extraction_summary.csv")
capture.output(sessionInfo(), file="metadata/outcome_coloc_extraction_R_sessionInfo.txt")
cat("PASS: all primary outcomes and strict-AS phenotype sensitivity extracted\n")
