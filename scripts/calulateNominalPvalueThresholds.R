suppressMessages(library(qvalue))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
  cat("Usage: Rscript calulateNominalPvalueThresholds.R <perm_file> <fdr> <output_file>\n\n")
  cat("Arguments:\n")
  cat("  perm_file    Permutation pass output from rust-fastqtl (plain or .gz)\n")
  cat("  fdr          FDR threshold (e.g. 0.05 for 5%)\n")
  cat("  output_file  Output file: phenotype_id, nominal_threshold, q-value\n\n")
  cat("Example:\n")
  cat("  Rscript calulateNominalPvalueThresholds.R perm.txt.gz 0.05 perm.fdr05.txt.gz\n")
  quit(status = 0)
}

if (length(args) < 3) {
  cat("Error: 3 arguments required. Run with -h for usage.\n")
  quit(status = 1)
}

ifile = args[1]
fdr = as.numeric(args[2]);

cat("Processing fastQTL concatenated output [", ifile, "] controlling for FDR =", fdr * 100, "%\n");

#Read data
D = read.table(ifile, hea=FALSE, stringsAsFactors=FALSE)
# Filter: remove rows with missing/infinite/out-of-range beta-adjusted p-values
D = D[which(is.finite(D[, 11]) & D[, 11] > 0 & D[, 11] <= 1), ]
cat("  * Number of molecular phenotypes =" , nrow(D), "\n")
cat("  * Correlation between Beta approx. and Empirical p-values =", round(cor(D[, 9], D[, 10]), 4), "\n")

#Run qvalue on pvalues for best signals
Q = tryCatch(
  qvalue(D[, 11]),
  error = function(e1) tryCatch(
    { message("  * smooth.spline failed, trying restricted lambda range")
      qvalue(D[, 11], lambda = seq(0.05, 0.85, 0.05)) },
    error = function(e2) tryCatch(
      { message("  * restricted lambda failed, trying bootstrap method")
        qvalue(D[, 11], pi0.method = "bootstrap") },
      error = function(e3) {
        message("  * all pi0 estimation methods failed, using pi0=1 (conservative BH-equivalent)")
        qvalue(D[, 11], pi0 = 1)
      }
    )
  )
)
D$qval = Q$qvalue
cat("  * Proportion of significant phenotypes =" , round((1 - Q$pi0) * 100, 2), "%\n")

#Determine significance threshold
set0 = D[which(D$qval <= fdr),] 
set1 = D[which(D$qval > fdr),]
pthreshold = (sort(set1$V11)[1] - sort(-1.0 * set0$V11)[1]) / 2
cat("  * Corrected p-value threshold = ", pthreshold, "\n")

#Calculate nominal pvalue thresholds
D$nthresholds = qbeta(pthreshold, D$V3, D$V4, ncp = 0, lower.tail = TRUE, log.p = FALSE)

#Write output
out_con <- if (grepl("\\.gz$", args[3])) gzfile(args[3], "w") else args[3]
write.table(D[, c(1, 13, 12)], out_con, quote=FALSE, row.names=FALSE, col.names=FALSE)
if (inherits(out_con, "connection")) close(out_con)

cat("Done\n")
