options(warn=1)
stopifnot(requireNamespace("Rsamtools",quietly=TRUE),requireNamespace("GenomicRanges",quietly=TRUE))
project <- "C:/Users/Lenovo/IBD_rebuild_runtime/source_project"
runtime <- "C:/Users/Lenovo/IBD_rebuild_runtime/finngen_R12"
sources <- c(
 M13_ANKYLOSPON=file.path(project,"data_external/finngen_as_r12/finngen_R12_M13_ANKYLOSPON.gz"),
 M13_ANKYLOSPON_STRICT=file.path(project,"data_external/finngen_as_r12_strict/finngen_R12_M13_ANKYLOSPON_STRICT.gz"),
 H7_IRIDOCYCLITIS=file.path(runtime,"finngen_R12_H7_IRIDOCYCLITIS.gz"),
 L12_PSORIASIS=file.path(runtime,"finngen_R12_L12_PSORIASIS.gz"),
 M13_PSORIARTH=file.path(runtime,"finngen_R12_M13_PSORIARTH.gz")
)
coords <- read.csv("tables/instrument_coordinate_audit.csv",stringsAsFactors=FALSE)
stopifnot(nrow(coords)==218L,all(coords$reciprocal_exact),!anyDuplicated(coords$MarkerName))
ranges <- GenomicRanges::GRanges(sub("^chr","",coords$chrom38),IRanges::IRanges(coords$pos38,width=1),strand="*")
header <- c("chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
result <- list(); k <- 0L
for(endpoint in names(sources)) {
  path <- sources[[endpoint]]
  stopifnot(file.exists(path),file.exists(paste0(path,".tbi")))
  lines <- Rsamtools::scanTabix(Rsamtools::TabixFile(path),param=ranges)
  stopifnot(length(lines)==nrow(coords))
  for(i in seq_along(lines)) for(line in lines[[i]]) {
    values <- strsplit(line,"\t",fixed=TRUE)[[1]]
    stopifnot(length(values)==length(header))
    k <- k+1L
    result[[k]] <- c(endpoint=endpoint,query_marker=coords$MarkerName[i],query_rsid=coords$reference_rsid[i],
                     query_chrom38=coords$chrom38[i],query_pos38=as.character(coords$pos38[i]),setNames(values,header))
  }
  cat(endpoint,"queries",length(lines),"returned rows",sum(lengths(lines)),"\n")
}
raw <- as.data.frame(do.call(rbind,result),stringsAsFactors=FALSE)
write.csv(raw,"tables/outcome_instrument_point_query_raw.csv",row.names=FALSE,na="")
capture.output(sessionInfo(),file="metadata/outcome_point_extraction_R_sessionInfo.txt")
cat("PASS: point query rows",nrow(raw),"\n")
