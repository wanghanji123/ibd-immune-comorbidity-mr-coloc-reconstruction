options(warn=1)
ns <- loadNamespace("TwoSampleMR")
fixtures <- list(
  underdispersed=list(bx=c(.1,.2,.3,.4), by=c(.0301,.0599,.0901,.1199), sy=rep(.05,4)),
  heterogeneous=list(bx=c(.1,.2,.3,.4), by=c(.3,-.1,.4,.05), sy=c(.02,.03,.04,.02)),
  two_snps=list(bx=c(.1,.3), by=c(.03,.05), sy=c(.03,.04)),
  signed=list(bx=c(-.1,.2,-.3,.4), by=c(-.0301,.0599,-.0901,.1199), sy=rep(.05,4))
)
out <- do.call(rbind, lapply(names(fixtures), function(n) {
  d <- fixtures[[n]]
  fit <- get("mr_ivw",ns)(d$bx,d$by,rep(.01,length(d$bx)),d$sy)
  data.frame(fixture=n,b=fit$b,se=fit$se,pval=fit$pval,Q=fit$Q,Q_df=fit$Q_df,Q_pval=fit$Q_pval)
}))
write.csv(out,"tables/synthetic_MR_R_crosscheck.csv",row.names=FALSE)
w <- get("mr_wald_ratio",ns)(.2,.1,.01,.02,get("default_parameters",ns)())
stopifnot(abs(w$b-.5)<1e-12,abs(w$se-.1)<1e-12,w$nsnp==1)
adj <- p.adjust(c(.01,rep(1,11)),method="BH")
stopifnot(abs(adj[1]-.12)<1e-12)
capture.output(sessionInfo(),file="metadata/MR_R_sessionInfo.txt")
cat("PASS: synthetic IVW, single-SNP Wald dispatch, fixed-family FDR; no real outcome MR run\n")
