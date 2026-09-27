pkgs <- c(
  "TwoSampleMR", "MRPRESSO", "coloc", "VariantAnnotation",
  "rtracklayer", "data.table"
)
installed <- installed.packages()
versions <- vapply(
  pkgs,
  function(pkg) {
    if (pkg %in% rownames(installed)) as.character(installed[pkg, "Version"]) else NA_character_
  },
  character(1)
)
out <- data.frame(package = pkgs, version = versions)
write.csv(out, "environment/R_key_package_versions.csv", row.names = FALSE)
