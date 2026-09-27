options(warn=1)
ns <- loadNamespace("TwoSampleMR")
param <- get("default_parameters",ns)(); param$nboot <- 10000L
exposures <- c(GCST004131="IBD",GCST004132="CD",GCST004133="UC")
branches <- c("primary_all","primary_nonMHC","relaxed_cohort_p_all","relaxed_cohort_p_nonMHC")
outcomes <- c("M13_ANKYLOSPON","H7_IRIDOCYCLITIS","L12_PSORIASIS","M13_PSORIARTH","M13_ANKYLOSPON_STRICT")
primary.outcomes <- outcomes[1:4]
results <- list(); loo <- list(); k <- 0L; j <- 0L
to.row <- function(fit,method,exposure,id,branch,outcome,seed) {
 data.frame(exposure=exposure,gwas_id=id,branch=branch,outcome=outcome,method=method,seed=seed,
   nsnp=fit$nsnp,beta=fit$b,se=fit$se,p=fit$pval,OR=exp(fit$b),CI_low=exp(fit$b-1.96*fit$se),CI_high=exp(fit$b+1.96*fit$se),
   Q=if(is.null(fit$Q)) NA else fit$Q,Q_df=if(is.null(fit$Q_df)) NA else fit$Q_df,
   Q_p=if(is.null(fit$Q_pval)) NA else fit$Q_pval,
   egger_intercept=if(is.null(fit$b_i)) NA else fit$b_i,
   egger_intercept_se=if(is.null(fit$se_i)) NA else fit$se_i,
   egger_intercept_p=if(is.null(fit$pval_i)) NA else fit$pval_i,stringsAsFactors=FALSE)
}
for(id in names(exposures)) for(branch in branches) for(outcome in outcomes) {
 path <- file.path("tables",paste0("harmonised_",id,"_",branch,"_to_",outcome,".csv"))
 d <- read.csv(path,stringsAsFactors=FALSE)
 stopifnot(nrow(d)>=2L,!anyDuplicated(d$MarkerName),all(is.finite(d$beta_exposure)),all(is.finite(d$se_exposure)),
           all(is.finite(d$beta_outcome)),all(is.finite(d$se_outcome)),all(d$se_exposure>0),all(d$se_outcome>0))
 bx<-d$beta_exposure; by<-d$beta_outcome; sx<-d$se_exposure; sy<-d$se_outcome
 seed <- 270000L + match(id,names(exposures))*1000L + match(branch,branches)*100L + match(outcome,outcomes)
 set.seed(seed); ivw <- get("mr_ivw",ns)(bx,by,sx,sy,param)
 set.seed(seed); med <- get("mr_weighted_median",ns)(bx,by,sx,sy,param)
 set.seed(seed); egger <- get("mr_egger_regression",ns)(bx,by,sx,sy,param)
 k<-k+1L; results[[k]]<-to.row(ivw,"IVW multiplicative random effects",exposures[[id]],id,branch,outcome,seed)
 k<-k+1L; results[[k]]<-to.row(med,"Weighted median",exposures[[id]],id,branch,outcome,seed)
 k<-k+1L; results[[k]]<-to.row(egger,"MR-Egger",exposures[[id]],id,branch,outcome,seed)
 if(branch=="primary_all" && outcome %in% primary.outcomes) for(i in seq_len(nrow(d))) {
   keep <- seq_len(nrow(d)) != i
   fit <- get("mr_ivw",ns)(bx[keep],by[keep],sx[keep],sy[keep],param)
   j<-j+1L; loo[[j]]<-data.frame(exposure=exposures[[id]],gwas_id=id,outcome=outcome,omitted=d$MarkerName[i],
     nsnp=fit$nsnp,beta=fit$b,se=fit$se,p=fit$pval,OR=exp(fit$b),CI_low=exp(fit$b-1.96*fit$se),CI_high=exp(fit$b+1.96*fit$se))
 }
}
res<-do.call(rbind,results)
primary.idx <- res$branch=="primary_all" & res$outcome %in% primary.outcomes & res$method=="IVW multiplicative random effects"
stopifnot(sum(primary.idx)==12L)
res$FDR_BH_primary12 <- NA_real_
res$FDR_BH_primary12[primary.idx] <- p.adjust(res$p[primary.idx],method="BH",n=12L)
write.csv(res,"tables/mr_all_locked_branches.csv",row.names=FALSE,na="")
write.csv(res[primary.idx,],"tables/mr_primary12.csv",row.names=FALSE,na="")
write.csv(do.call(rbind,loo),"tables/mr_primary12_leave_one_out.csv",row.names=FALSE,na="")
capture.output(sessionInfo(),file="metadata/actual_MR_R_sessionInfo.txt")
cat("PASS: MR estimates",nrow(res),"primary IVW hypotheses",sum(primary.idx),"leave-one-out rows",length(loo),"\n")
