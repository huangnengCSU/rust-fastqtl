#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(qvalue))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript calculateNominalPvalueThresholds.R <perm_file> <fdr> <output_file>\n")
  quit(status = 1)
}

ifile <- args[1]
fdr   <- as.numeric(args[2])
ofile <- args[3]

cat("Processing permutation output:", ifile, "\n")
cat("Target FDR:", fdr, "\n")

D <- read.table(ifile, header = FALSE, stringsAsFactors = FALSE)

if (ncol(D) != 17) {
  stop("Expected 17 columns in permutation file, found ", ncol(D))
}

colnames(D) <- c(
  "gene_id", "num_var", "beta_shape1", "beta_shape2", "true_df",
  "pval_true_df", "variant_id", "tss_distance", "ma_samples", "ma_count",
  "maf", "ref_factor", "pval_nominal", "slope", "slope_se",
  "pval_perm", "pval_beta"
)

# keep valid pval_beta
D <- D[is.finite(D$pval_beta) & D$pval_beta > 0 & D$pval_beta <= 1, ]

cat("  * Number of genes =", nrow(D), "\n")
cat("  * Correlation between beta and empirical p-values =",
    round(cor(D$pval_beta, D$pval_perm), 4), "\n")

Q <- tryCatch(
  qvalue(D$pval_beta),
  error = function(e1) tryCatch(
    {
      message("  * smooth.spline failed, trying restricted lambda range")
      qvalue(D$pval_beta, lambda = seq(0.05, 0.85, 0.05))
    },
    error = function(e2) tryCatch(
      {
        message("  * restricted lambda failed, trying bootstrap")
        qvalue(D$pval_beta, pi0.method = "bootstrap")
      },
      error = function(e3) {
        message("  * all methods failed, using pi0=1")
        qvalue(D$pval_beta, pi0 = 1)
      }
    )
  )
)

D$qval <- Q$qvalues

sig <- D[D$qval <= fdr, ]
nonsig <- D[D$qval > fdr, ]

if (nrow(sig) == 0) {
  stop("No genes pass the chosen FDR.")
}

if (nrow(nonsig) == 0) {
  pthreshold <- max(sig$pval_beta, na.rm = TRUE)
} else {
  lb <- max(sig$pval_beta, na.rm = TRUE)
  ub <- min(nonsig$pval_beta, na.rm = TRUE)
  pthreshold <- (lb + ub) / 2
}

cat("  * Corrected p-value threshold =", pthreshold, "\n")

D$pval_nominal_threshold <- qbeta(
  pthreshold,
  shape1 = D$beta_shape1,
  shape2 = D$beta_shape2
)

out <- D[, c("gene_id", "pval_nominal_threshold", "qval")]

out_con <- if (grepl("\\.gz$", ofile)) gzfile(ofile, "w") else ofile
write.table(out, out_con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
if (inherits(out_con, "connection")) close(out_con)

cat("Done\n")