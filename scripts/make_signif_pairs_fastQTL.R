#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript make_signif_pairs.R <threshold_file> <nominal_file> <output_file>\n")
  quit(status = 1)
}

threshold_file <- args[1]
nominal_file   <- args[2]
out_file       <- args[3]

message("[1/4] Reading threshold file: ", threshold_file)
thr <- fread(threshold_file)

required_thr <- c("gene_id", "pval_nominal_threshold", "qval")
miss_thr <- setdiff(required_thr, names(thr))
if (length(miss_thr) > 0) {
  stop("Missing columns in threshold file: ", paste(miss_thr, collapse = ", "))
}

thr_sig <- thr[qval <= 0.05]

message("[2/4] Reading nominal file: ", nominal_file)
nom <- fread(nominal_file, header = FALSE)

if (ncol(nom) != 9) {
  stop("Expected 9 columns in nominal file, found ", ncol(nom))
}

setnames(nom, c(
  "gene_id", "variant_id", "tss_distance", "ma_samples",
  "ma_count", "maf", "pval_nominal", "slope", "slope_se"
))

message("[3/4] Merging and filtering")
res <- merge(nom, thr_sig, by = "gene_id", all = FALSE)
sig_pairs <- res[pval_nominal <= pval_nominal_threshold]

setcolorder(sig_pairs, c(
  "gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count",
  "maf", "pval_nominal", "pval_nominal_threshold", "qval",
  "slope", "slope_se"
))

message("[4/4] Writing output: ", out_file)
fwrite(sig_pairs, out_file, sep = "\t", quote = FALSE)

message("Done. Significant pairs: ", nrow(sig_pairs))