options(warn=1)
stopifnot(requireNamespace("data.table", quietly=TRUE))
stopifnot(requireNamespace("rtracklayer", quietly=TRUE))
stopifnot(requireNamespace("GenomicRanges", quietly=TRUE))
stopifnot(requireNamespace("IRanges", quietly=TRUE))

input <- "C:/Users/Lenovo/IBD_rebuild_runtime/round03_coloc/exposure_regions_GRCh37.csv.gz"
runtime <- "C:/Users/Lenovo/IBD_rebuild_runtime/round03_coloc"
dir.create(runtime, recursive=TRUE, showWarnings=FALSE)
forward <- rtracklayer::import.chain("metadata/hg19ToHg38.over.chain")
reverse <- rtracklayer::import.chain("metadata/hg38ToHg19.over.chain")
dat <- data.table::fread(input, showProgress=TRUE)
dat[, source_row := .I]
dat[, `:=`(forward_mappings=0L, chrom38=NA_character_, pos38=NA_integer_,
           mapping_strand=NA_character_, reciprocal_exact=FALSE)]

chunk_size <- 100000L
for (gwas in unique(dat$gwas_id)) {
  ids <- which(dat$gwas_id == gwas)
  chunks <- split(ids, ceiling(seq_along(ids) / chunk_size))
  for (ids_chunk in chunks) {
    gr <- GenomicRanges::GRanges(paste0("chr", dat$chrom37[ids_chunk]),
      IRanges::IRanges(dat$pos37[ids_chunk], width=1), strand="+")
    lifted <- rtracklayer::liftOver(gr, forward)
    counts <- lengths(lifted)
    dat$forward_mappings[ids_chunk] <- counts
    one_local <- which(counts == 1L)
    if (!length(one_local)) next
    mapped <- unlist(lifted[one_local], use.names=FALSE)
    one_global <- ids_chunk[one_local]
    back <- rtracklayer::liftOver(mapped, reverse)
    back_count <- lengths(back)
    back_one_local <- which(back_count == 1L)
    reciprocal <- rep(FALSE, length(one_local))
    if (length(back_one_local)) {
      back_gr <- unlist(back[back_one_local], use.names=FALSE)
      original <- gr[one_local[back_one_local]]
      reciprocal[back_one_local] <-
        as.character(GenomicRanges::seqnames(back_gr)) == as.character(GenomicRanges::seqnames(original)) &
        IRanges::start(back_gr) == IRanges::start(original) &
        IRanges::width(back_gr) == 1L &
        as.character(GenomicRanges::strand(back_gr)) == "+"
    }
    dat$chrom38[one_global] <- sub("^chr", "", as.character(GenomicRanges::seqnames(mapped)))
    dat$pos38[one_global] <- IRanges::start(mapped)
    dat$mapping_strand[one_global] <- as.character(GenomicRanges::strand(mapped))
    dat$reciprocal_exact[one_global] <- reciprocal
    cat(gwas, "chunk", min(ids_chunk), max(ids_chunk), "single", length(one_local),
        "reciprocal", sum(reciprocal), "\n")
  }
}

dat[, ref38 := ref37]
dat[, alt38 := alt37]
dat[, effect_allele38 := effect_allele37]
dat[, other_allele38 := other_allele37]
rev <- which(dat$mapping_strand == "-")
for (field in c("ref38", "alt38", "effect_allele38", "other_allele38"))
  data.table::set(dat, i=rev, j=field, value=chartr("ACGT", "TGCA", dat[[field]][rev]))
dat[, variant_key := ifelse(reciprocal_exact,
  paste(chrom38, pos38, ref38, alt38, sep=":"), NA_character_)]
dat[, duplicate_mapped_key := reciprocal_exact &
  (duplicated(paste(gwas_id, region_id, variant_key, sep="|")) |
   duplicated(paste(gwas_id, region_id, variant_key, sep="|"), fromLast=TRUE))]
dat[, mapping_status := data.table::fcase(
  forward_mappings != 1L, "HOLD_NONUNIQUE_FORWARD_MAPPING",
  !reciprocal_exact, "HOLD_RECIPROCAL_MAPPING_FAILURE",
  duplicate_mapped_key, "HOLD_DUPLICATE_MAPPED_KEY",
  default="PASS_RECIPROCAL_UNIQUE_MAPPING")]

audit_path <- file.path(runtime, "exposure_coloc_liftover_audit.csv.gz")
pass_path <- file.path(runtime, "exposure_coloc_variants_GRCh38.csv.gz")
data.table::fwrite(dat, audit_path, compress="gzip")
pass <- dat[mapping_status == "PASS_RECIPROCAL_UNIQUE_MAPPING"]
stopifnot(!anyDuplicated(pass[, .(gwas_id, region_id, variant_key)]))
data.table::fwrite(pass, pass_path, compress="gzip")
summary <- dat[, .(input_variants=.N, reciprocal_unique_pass=sum(mapping_status == "PASS_RECIPROCAL_UNIQUE_MAPPING"),
  nonunique_forward=sum(mapping_status == "HOLD_NONUNIQUE_FORWARD_MAPPING"),
  reciprocal_failure=sum(mapping_status == "HOLD_RECIPROCAL_MAPPING_FAILURE"),
  duplicate_mapped_key=sum(mapping_status == "HOLD_DUPLICATE_MAPPED_KEY"),
  reverse_strand=sum(mapping_strand == "-", na.rm=TRUE)), by=.(trait, gwas_id)]
data.table::fwrite(summary, "tables/exposure_coloc_liftover_summary.csv")
capture.output(sessionInfo(), file="metadata/coloc_liftover_R_sessionInfo.txt")
cat("Wrote", audit_path, "and", pass_path, "\n")
