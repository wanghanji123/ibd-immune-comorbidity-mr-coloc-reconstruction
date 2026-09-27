options(warn=1)
stopifnot(requireNamespace("rtracklayer", quietly=TRUE))
forward <- rtracklayer::import.chain("metadata/hg19ToHg38.over.chain")
reverse <- rtracklayer::import.chain("metadata/hg38ToHg19.over.chain")
dat <- read.csv("tables/unique_instrument_variants_GRCh37.csv", stringsAsFactors=FALSE)
gr <- GenomicRanges::GRanges(paste0("chr", dat$chrom37), IRanges::IRanges(dat$pos37, width=1), strand="+")
lift <- rtracklayer::liftOver(gr, forward)
dat$forward_mappings <- lengths(lift)
dat$chrom38 <- NA_character_
dat$pos38 <- NA_integer_
dat$mapping_strand <- NA_character_
dat$reciprocal_exact <- FALSE
for (i in which(lengths(lift) == 1L)) {
  point <- lift[[i]]
  if (IRanges::width(point) != 1L) next
  back <- rtracklayer::liftOver(point, reverse)[[1]]
  dat$chrom38[i] <- as.character(GenomicRanges::seqnames(point))
  dat$pos38[i] <- IRanges::start(point)
  dat$mapping_strand[i] <- as.character(GenomicRanges::strand(point))
  dat$reciprocal_exact[i] <- length(back) == 1L &&
    as.character(GenomicRanges::seqnames(back)) == as.character(GenomicRanges::seqnames(gr[i])) &&
    IRanges::start(back) == IRanges::start(gr[i]) && IRanges::width(back) == 1L &&
    as.character(GenomicRanges::strand(back)) == "+"
}
dat$effect_allele38 <- dat$effect_allele
dat$other_allele38 <- dat$other_allele
reversed <- which(dat$mapping_strand == "-")
dat$effect_allele38[reversed] <- chartr("ACGT", "TGCA", dat$effect_allele[reversed])
dat$other_allele38[reversed] <- chartr("ACGT", "TGCA", dat$other_allele[reversed])
dat$mapping_status <- ifelse(dat$reciprocal_exact, "RECIPROCAL_COORDINATE_PASS_OUTCOME_ALLELES_NOT_YET_CHECKED", "HOLD_MAPPING")
write.csv(dat, "tables/instrument_coordinate_audit.csv", row.names=FALSE, na="")
cat("variants:", nrow(dat), " reciprocal pass:", sum(dat$reciprocal_exact), " reverse strand:", length(reversed), "\n")
capture.output(sessionInfo(), file="metadata/liftover_R_sessionInfo.txt")
