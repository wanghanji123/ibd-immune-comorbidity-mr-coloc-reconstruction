args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L)
id<-args[1]; exposure<-args[2]; outcome<-args[3]; seed<-as.integer(args[4]); n_sim<-as.integer(args[5])
input<-file.path("tables",paste0("harmonised_",id,"_primary_all_to_",outcome,".csv"))
outdir<-file.path("results","mrpresso",paste(id,outcome,paste0("seed",seed),paste0("n",n_sim),sep="__"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(input,stringsAsFactors=FALSE)
stopifnot(nrow(d)>=2L,!anyDuplicated(d$MarkerName),n_sim>nrow(d)/.05)
warnings_seen<-character()
start<-Sys.time()
fit<-withCallingHandlers(
  MRPRESSO::mr_presso(BetaOutcome="beta_outcome",BetaExposure="beta_exposure",SdOutcome="se_outcome",SdExposure="se_exposure",
    data=d,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,SignifThreshold=.05,NbDistribution=n_sim,seed=seed),
  warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
elapsed<-as.numeric(difftime(Sys.time(),start,units="secs"))
saveRDS(fit,file.path(outdir,"raw_result.rds"))
main<-fit[["Main MR results"]]
main$exposure<-exposure; main$gwas_id<-id; main$outcome<-outcome; main$seed<-seed; main$n_sim<-n_sim; main$nsnp_input<-nrow(d)
write.csv(main,file.path(outdir,"main_results.csv"),row.names=FALSE,na="")
global<-fit[["MR-PRESSO results"]][["Global Test"]]
write.csv(data.frame(exposure=exposure,gwas_id=id,outcome=outcome,seed=seed,n_sim=n_sim,nsnp_input=nrow(d),
  RSSobs=global$RSSobs,Pvalue=as.character(global$Pvalue),elapsed_seconds=elapsed),file.path(outdir,"global_test.csv"),row.names=FALSE)
outlier<-fit[["MR-PRESSO results"]][["Outlier Test"]]
if(!is.null(outlier)) {
  stopifnot(nrow(outlier)==nrow(d)); outlier$MarkerName<-d$MarkerName
  write.csv(outlier,file.path(outdir,"outlier_test.csv"),row.names=FALSE,na="")
}
distortion<-fit[["MR-PRESSO results"]][["Distortion Test"]]
if(!is.null(distortion)) {
  indices<-distortion[["Outliers Indices"]]
  write.csv(data.frame(exposure=exposure,gwas_id=id,outcome=outcome,seed=seed,n_sim=n_sim,
    outlier_indices=paste(indices,collapse=";"),outlier_markers=paste(d$MarkerName[indices],collapse=";"),
    distortion_coefficient=as.numeric(distortion[["Distortion Coefficient"]]),
    distortion_p=as.character(distortion[["Pvalue"]])),file.path(outdir,"distortion_test.csv"),row.names=FALSE,na="")
}
writeLines(unique(warnings_seen),file.path(outdir,"warnings.txt"))
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
cat(outdir,"elapsed",elapsed,"seconds warnings",length(unique(warnings_seen)),"\n")
