#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript calulateNominalPvalueThresholds_rustQTL.R <perm_file> <fdr> <output_file>\n")
  quit(status = 1)
}

ifile <- args[1]
fdr <- as.numeric(args[2])
ofile <- args[3]

read_perm_block_format <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) %% 3 != 0) {
    stop("Permutation block-format file should have 3 lines per phenotype; got ", length(lines), " non-empty lines")
  }

  split_ws <- function(x) strsplit(trimws(x), "\\s+")[[1]]

  parsed <- vector("list", length(lines) / 3)
  idx <- 1L

  for (i in seq(1L, length(lines), by = 3L)) {
    top <- split_ws(lines[i])
    mid <- split_ws(lines[i + 1L])
    bot <- split_ws(lines[i + 2L])

    if (length(top) != 6) {
      stop("Invalid top line in block ", idx, ": expected 6 fields, got ", length(top))
    }
    if (length(mid) != 9) {
      stop("Invalid middle line in block ", idx, ": expected 9 fields, got ", length(mid))
    }
    if (length(bot) != 2) {
      stop("Invalid bottom line in block ", idx, ": expected 2 fields, got ", length(bot))
    }

    parsed[[idx]] <- data.table(
      phenotype_id = top[1],
      n_variants = as.integer(top[2]),
      beta_shape1 = as.numeric(top[3]),
      beta_shape2 = as.numeric(top[4]),
      eff_df = as.numeric(top[5]),
      p_nom_eff = as.numeric(top[6]),
      variant_id = mid[1],
      distance = as.numeric(mid[2]),
      distance_to_body = as.numeric(mid[3]),
      ma_samples = as.integer(mid[4]),
      ma_count = as.integer(mid[5]),
      maf = as.numeric(mid[6]),
      p_nom = as.numeric(mid[7]),
      beta = as.numeric(mid[8]),
      beta_se = as.numeric(mid[9]),
      p_empirical = as.numeric(bot[1]),
      p_beta_adjusted = as.numeric(bot[2])
    )
    idx <- idx + 1L
  }

  rbindlist(parsed, use.names = TRUE)
}

read_perm_rust <- function(path) {
  tab <- suppressWarnings(fread(path, header = FALSE, sep = "\t", showProgress = FALSE))

  if (ncol(tab) == 17) {
    setnames(tab, c(
      "phenotype_id", "n_variants", "beta_shape1", "beta_shape2", "eff_df", "p_nom_eff",
      "variant_id", "distance", "distance_to_body", "ma_samples", "ma_count", "maf",
      "p_nom", "beta", "beta_se", "p_empirical", "p_beta_adjusted"
    ))
    return(tab)
  }

  if (ncol(tab) == 16) {
    setnames(tab, c(
      "phenotype_id", "n_variants", "beta_shape1", "beta_shape2", "eff_df", "p_nom_eff",
      "variant_id", "distance", "ma_samples", "ma_count", "maf",
      "p_nom", "beta", "beta_se", "p_empirical", "p_beta_adjusted"
    ))
    tab[, distance_to_body := NA_real_]
    setcolorder(tab, c(
      "phenotype_id", "n_variants", "beta_shape1", "beta_shape2", "eff_df", "p_nom_eff",
      "variant_id", "distance", "distance_to_body", "ma_samples", "ma_count", "maf",
      "p_nom", "beta", "beta_se", "p_empirical", "p_beta_adjusted"
    ))
    return(tab)
  }

  read_perm_block_format(path)
}

cat("Processing rust-fastqtl permutation output:", ifile, "\n")
cat("Target FDR:", fdr, "\n")

D <- read_perm_rust(ifile)

D <- D[is.finite(p_beta_adjusted) & p_beta_adjusted > 0 & p_beta_adjusted <= 1]

cat("  * Number of phenotypes =", nrow(D), "\n")
cat("  * Correlation between beta-adjusted and empirical p-values =",
    round(cor(D$p_beta_adjusted, D$p_empirical), 4), "\n")

Q <- tryCatch(
  qvalue(D$p_beta_adjusted),
  error = function(e1) tryCatch(
    {
      message("  * smooth.spline failed, trying restricted lambda range")
      qvalue(D$p_beta_adjusted, lambda = seq(0.05, 0.85, 0.05))
    },
    error = function(e2) tryCatch(
      {
        message("  * restricted lambda failed, trying bootstrap")
        qvalue(D$p_beta_adjusted, pi0.method = "bootstrap")
      },
      error = function(e3) {
        message("  * all methods failed, using pi0=1")
        qvalue(D$p_beta_adjusted, pi0 = 1)
      }
    )
  )
)

D[, qval := Q$qvalues]

sig <- D[qval <= fdr]
nonsig <- D[qval > fdr]

if (nrow(sig) == 0) {
  stop("No phenotypes pass the chosen FDR.")
}

if (nrow(nonsig) == 0) {
  pthreshold <- max(sig$p_beta_adjusted, na.rm = TRUE)
} else {
  lb <- max(sig$p_beta_adjusted, na.rm = TRUE)
  ub <- min(nonsig$p_beta_adjusted, na.rm = TRUE)
  pthreshold <- (lb + ub) / 2
}

cat("  * Corrected p-value threshold =", pthreshold, "\n")

D[, pval_nominal_threshold := qbeta(
  pthreshold,
  shape1 = beta_shape1,
  shape2 = beta_shape2
)]

out <- D[, .(phenotype_id, pval_nominal_threshold, qval)]

fwrite(out, ofile, sep = "\t", quote = FALSE)

cat("Done\n")
