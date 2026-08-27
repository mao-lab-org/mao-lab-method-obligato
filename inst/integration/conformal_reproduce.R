#!/usr/bin/env Rscript
################################################################################
# CONFORMAL REPRODUCTION GATE (S4 / Item 3).
#
# Confirms the packaged split-conformal functions (all_pair_ll,
# conformal_calibrate, conformal_sets) reproduce the reported 90%-coverage mean
# prediction-set sizes from Phase0_Detection/conformal_composition.R:
#     B2 ~6.7 | B9 ~7.6 | PBMC ~2.7.
#
# It replays that script's stratified 50/50 calibration/evaluation split (seed 42,
# ground-truth doublets) and its all-true-doublet nonconformity, but computes the
# log-likelihoods and the conformal set through the PACKAGE functions. This is the
# "results reproduce from the package alone" gate for conformal, analogous to M2b.
#
#   Rscript Obligato/inst/integration/conformal_reproduce.R batch9 full 7.6
################################################################################

.a    <- commandArgs(trailingOnly = TRUE)
LAB   <- if (length(.a) >= 1) .a[[1]] else "batch9"
RTAG  <- if (length(.a) >= 2) .a[[2]] else "full"
TARGET<- if (length(.a) >= 3) as.numeric(.a[[3]]) else NA_real_
SEED  <- 42L
ALPHA <- 0.10

# Run from the repo root (like the M2b gates); repo root is the working directory.
PROJ <- getwd()
stopifnot(dir.exists(file.path(PROJ, "Phase0_Detection")), dir.exists(file.path(PROJ, "Obligato")))
pkgload::load_all(file.path(PROJ, "Obligato"), quiet = TRUE)

.suffix  <- if (RTAG == "main") "" else paste0("_", RTAG)
rds_path <- file.path(PROJ, sprintf("Phase0_Detection/results/methodC_compose%s_%s.rds", .suffix, LAB))
es_path  <- file.path(PROJ, sprintf("Phase0_Detection/data/eval_set_%s.rds", LAB))

obj  <- readRDS(rds_path)
es   <- readRDS(es_path)
meta <- as.data.frame(es$test_metadata, stringsAsFactors = FALSE)
pc   <- obj$per_cell
stopifnot(nrow(meta) == nrow(pc))

# ---- all-pairs log-likelihoods via the PACKAGE function ----------------------
log_h <- all_pair_ll(obj$S_test, obj$pair_models)   # n x P, colnames = pair labels

# ---- ground-truth doublets + stratified 50/50 split (mirrors eval script) ----
dbl_idx <- which(pc$is_doublet == 1)
gt_pair <- pair_label(as.character(meta$type_A[dbl_idx]),
                      as.character(meta$type_B[dbl_idx]))

set.seed(SEED)
cal_local <- logical(length(dbl_idx))
for (pl in unique(gt_pair)) {
  ii <- which(gt_pair == pl)
  cal_local[sample(ii, floor(length(ii) / 2))] <- TRUE
}
cal_rows  <- dbl_idx[cal_local]
eval_rows <- dbl_idx[!cal_local]

# ---- calibrate + build sets via the PACKAGE functions ------------------------
q_hat <- conformal_calibrate(log_h[cal_rows, , drop = FALSE],
                             gt_pair[cal_local], alpha = ALPHA)[[1]]
cs    <- conformal_sets(log_h[eval_rows, , drop = FALSE], q_hat)

# coverage on the eval split
covered <- vapply(seq_along(eval_rows), function(i)
  gt_pair[!cal_local][i] %in% cs$sets[[i]], logical(1))

mean_sz <- mean(cs$sizes)
cat(sprintf("\n===== Conformal reproduction — %s [%s] =====\n", LAB, RTAG))
cat(sprintf("  q_hat(alpha=%.2f) = %.3f\n", ALPHA, q_hat))
cat(sprintf("  eval doublets     = %d\n", length(eval_rows)))
cat(sprintf("  empirical coverage= %.3f (target %.2f)\n", mean(covered), 1 - ALPHA))
cat(sprintf("  MEAN set size     = %.2f  (median %.0f)\n", mean_sz, median(cs$sizes)))
if (!is.na(TARGET)) {
  ok <- abs(mean_sz - TARGET) <= 0.5
  cat(sprintf("  target mean size  = %.2f  ->  %s (|diff| = %.2f, tol 0.5)\n",
              TARGET, if (ok) "PASS" else "FAIL", abs(mean_sz - TARGET)))
  if (!ok) quit(status = 1)
}
cat("Done.\n")
