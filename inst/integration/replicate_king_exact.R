################################################################################
# M2b PRIMARY GATE — exact reproduction of the King composition & detection from
# the STORED intermediates.
#
# The full-from-raw run (replicate_king_m2b.R) is confounded by RNG: PhiSpace/scran
# scoring is non-deterministic, and the doublet draw + seed differ from the stored
# run, so per-cell agreement across seeds is the wrong thing to demand. This gate
# instead isolates the package's COMPUTATION: given the stored scores and fitted
# models, do the package functions reproduce the stored per-cell outputs EXACTLY?
# No scoring, no training, no RNG — fast and deterministic.
################################################################################

suppressMessages({ library(xgboost); library(Matrix) })
PROJ <- Sys.getenv("OBLIGATO_PROJECT_ROOT", ".")
pkgload::load_all(file.path(PROJ, "Obligato"), quiet = TRUE)

o  <- readRDS(file.path(PROJ, "CaseStudies/results/kingetal_full_compose.rds"))
pc <- o$per_cell
S  <- o$S_all_cells

# --- 1. Composition: compose_pairs on stored scores + stored models -----------
comp <- compose_pairs(S, o$pair_models, candidates = "all5", output = "top1")
agree <- mapply(function(a, b, x, y) setequal(c(a, b), c(x, y)),
                comp$c1, comp$c2, pc$dbllik_all5_c1, pc$dbllik_all5_c2, USE.NAMES = FALSE)
comp_exact <- mean(agree)
# ordered exact (c1==c1 & c2==c2) too, in case tie-breaking is stricter
comp_ordered <- mean(comp$c1 == pc$dbllik_all5_c1 & comp$c2 == pc$dbllik_all5_c2)

# --- 2. Detection: compute_features + predict(stored classifier) ---------------
F_all <- compute_features(S, o$sing_models, o$pair_models)
det   <- predict(o$classifier, as.matrix(F_all))
det_maxdiff <- max(abs(det - pc$detection_score))
det_cor     <- cor(det, pc$detection_score)

cat("\n============ M2b KING — EXACT gate (stored intermediates) ============\n")
cat(sprintf("composition setequal-exact : %.6f  (n=%d)\n", comp_exact, nrow(comp)))
cat(sprintf("composition ordered-exact  : %.6f\n", comp_ordered))
cat(sprintf("detection  max|diff|       : %.3e\n", det_maxdiff))
cat(sprintf("detection  correlation     : %.8f\n", det_cor))

pass_comp <- comp_exact >= 0.9999
pass_det  <- det_maxdiff < 1e-5
cat(sprintf("\nCOMPOSITION EXACT: %s\nDETECTION EXACT:   %s\n",
            if (pass_comp) "PASS" else "FAIL",
            if (pass_det)  "PASS" else "FAIL"))

# Show a few disagreements if composition isn't exact, for diagnosis.
if (!pass_comp) {
  bad <- which(!agree)[1:min(8, sum(!agree))]
  cat("\nfirst composition disagreements (pkg vs stored):\n")
  print(data.frame(cell = comp$cell[bad],
                   pkg = paste(comp$c1[bad], comp$c2[bad], sep = "+"),
                   stored = paste(pc$dbllik_all5_c1[bad], pc$dbllik_all5_c2[bad], sep = "+")))
}

saveRDS(list(comp_exact = comp_exact, comp_ordered = comp_ordered,
             det_maxdiff = det_maxdiff, det_cor = det_cor,
             pass_comp = pass_comp, pass_det = pass_det),
        file.path(PROJ, "CaseStudies/results/m2b_king_exact.rds"))
cat("\nM2B_EXACT_DONE\n")
