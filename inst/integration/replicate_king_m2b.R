################################################################################
# M2b REPLICATION GATE — King et al. all-cells composition, package alone.
#
# Runs Obligato::compose() end-to-end on the King object and compares to the
# stored pipeline output (CaseStudies/results/kingetal_full_compose_percell.csv).
#
# Per the determinism decision (package_plan Phase B/C), reproduction of the OLD
# stored numbers is TOLERANCE-based, not exact: detection inherits xgboost RNG and
# composition inherits the synthetic-doublet sampling RNG. The gate therefore
# checks that a deterministic package run lands within tolerance of the stored
# draw. Re-running this script twice must give byte-identical package output
# (that exactness is asserted at the end).
#
# Not part of the fast testthat suite (needs PhiSpace + the 284 MB King object).
# Run:  OBLIGATO_PROJECT_ROOT=<proj> Rscript inst/integration/replicate_king_m2b.R
################################################################################

suppressMessages({
  library(SingleCellExperiment); library(PhiSpace); library(xgboost); library(Matrix)
})

PROJ <- Sys.getenv("OBLIGATO_PROJECT_ROOT", ".")
PKG  <- Sys.getenv("OBLIGATO_PKG", file.path(PROJ, "Obligato"))
pkgload::load_all(PKG, quiet = TRUE)

t0 <- Sys.time()
say <- function(...) cat(sprintf("[%5.0fs] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n")

# --- Acceptance tolerances (documented; adjust with justification) ------------
TOL <- list(
  composition_agreement = 0.95,   # frac cells with setequal all5 pair vs stored
  detection_spearman    = 0.98,   # rank agreement of detection score vs stored
  flagged_jaccard       = 0.90,   # overlap of flagged sets at each own threshold
  threshold_abs         = 0.05    # |pkg threshold - stored 0.711|
)

# --- Inputs -------------------------------------------------------------------
say("loading King object + B9 reference")
k   <- readRDS(file.path(PROJ, "CaseStudies/data/kingetal_sce_full.rds"))
ref <- readRDS(file.path(PROJ, "CaseStudies/data/ref_tonsil_b9.rds"))
stored <- read.csv(file.path(PROJ, "CaseStudies/results/kingetal_full_compose_percell.csv"),
                   stringsAsFactors = FALSE)

query     <- counts(k)
hc        <- as.logical(k$hc_singlet)
hc_labels <- as.character(k$phi_top1_azimuth)[hc]      # King's precomputed top-1
cell_id   <- colnames(k)
say(sprintf("query %d cells x %d genes | hc singlets %d | ref types %d",
            ncol(query), nrow(query), sum(hc), length(unique(ref[["predicted.celltype.l1"]]))))

# --- Run the package ----------------------------------------------------------
say("compose() run 1 (this scores the joint training pool via PhiSpace — slow)")
res <- compose(query, ref, phenotypes = "predicted.celltype.l1",
               hc_singlets = hc, hc_labels = hc_labels,
               detector = "builtin", n_per_pair = 200, seed = 1L)
say(sprintf("compose done: threshold=%.4f  flagged=%d (%.1f%%)",
            res$info$threshold, sum(res$flag), 100 * mean(res$flag)))

comp <- res$composition
comp$cell_id <- cell_id[comp$cell]
comp$det     <- res$detection_score
comp$flag    <- res$flag

# --- Align to stored ----------------------------------------------------------
m <- merge(comp, stored, by = "cell_id", suffixes = c("", "_stored"))
say(sprintf("aligned %d / %d cells to stored CSV", nrow(m), nrow(comp)))
stopifnot(nrow(m) > 0.99 * nrow(comp))

# --- Metrics ------------------------------------------------------------------
agree <- mapply(function(a, b, x, y) setequal(c(a, b), c(x, y)),
                m$c1, m$c2, m$dbllik_all5_c1, m$dbllik_all5_c2, USE.NAMES = FALSE)
comp_agreement <- mean(agree)

det_spearman <- suppressWarnings(cor(m$det, m$detection_score, method = "spearman"))
det_maxdiff  <- max(abs(m$det - m$detection_score))

stored_thr  <- 0.7100192   # regenerated 7-feature King reference (S2/Item 1, 2026-08-25)
pkg_flag    <- m$det > res$info$threshold
stored_flag <- m$detection_score > stored_thr
jac <- sum(pkg_flag & stored_flag) / sum(pkg_flag | stored_flag)
thr_diff <- abs(res$info$threshold - stored_thr)

# --- Verdict ------------------------------------------------------------------
checks <- data.frame(
  metric    = c("composition_agreement", "detection_spearman", "flagged_jaccard", "threshold_abs"),
  value     = c(comp_agreement, det_spearman, jac, thr_diff),
  tolerance = c(TOL$composition_agreement, TOL$detection_spearman, TOL$flagged_jaccard, TOL$threshold_abs),
  pass      = c(comp_agreement >= TOL$composition_agreement,
                det_spearman   >= TOL$detection_spearman,
                jac            >= TOL$flagged_jaccard,
                thr_diff       <= TOL$threshold_abs))
cat("\n==================== M2b KING RESULT ====================\n")
print(checks, row.names = FALSE)
cat(sprintf("detection max|diff| = %.4f (informational)\n", det_maxdiff))
overall <- all(checks$pass)
cat(sprintf("OVERALL: %s\n", if (overall) "PASS (within tolerance of stored draw)" else "FAIL"))

# --- Self-determinism: run 2 must be byte-identical ---------------------------
say("compose() run 2 (determinism check)")
res2 <- compose(query, ref, phenotypes = "predicted.celltype.l1",
                hc_singlets = hc, hc_labels = hc_labels,
                detector = "builtin", n_per_pair = 200, seed = 1L)
det_identical  <- identical(res$detection_score, res2$detection_score)
comp_identical <- identical(res$composition, res2$composition)
cat(sprintf("SELF-DETERMINISM: detection identical=%s  composition identical=%s\n",
            det_identical, comp_identical))

# --- Persist ------------------------------------------------------------------
out <- list(checks = checks, overall = overall, det_maxdiff = det_maxdiff,
            self_determinism = c(detection = det_identical, composition = comp_identical),
            pkg_threshold = res$info$threshold, n_aligned = nrow(m),
            timestamp = Sys.time())
saveRDS(out, file.path(PROJ, "CaseStudies/results/m2b_king_result.rds"))
write.csv(checks, file.path(PROJ, "CaseStudies/results/m2b_king_checks.csv"), row.names = FALSE)
say("saved m2b_king_result.rds + m2b_king_checks.csv")
cat("M2B_KING_DONE\n")
