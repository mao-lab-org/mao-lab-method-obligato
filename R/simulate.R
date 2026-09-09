##############################################################################
# Synthetic training doublets (M1b).
#
# Extracted from run_compose_kingetal.R:184-201 (Section 1), identical to
# run_methodC_compose.R's process_batch simulation block. For every heterotypic
# type pair, sample `n_per_pair` high-confidence singlets of each type and sum
# their raw counts into a synthetic doublet tagged with the unordered pair label.
#
# DETERMINISM (same sanctioned decision as train_detector, 2026-08-20): the scripts
# relied on the ambient RNG stream seeded once at the top; a package function
# cannot, so it seeds internally from `seed`. Consequence: the synthetic doublets —
# and therefore the fitted pair models and every composition call downstream — are
# reproducible run-to-run but do NOT bit-match the original stored draw. M2b's
# reproduction of the OLD stored numbers is therefore tolerance-based for
# composition as well as detection; the package reproduces ITSELF exactly.
#
# Behaviour-preserving deviations: N_PER_PAIR -> `n_per_pair` (locked default 200);
# ambient `types` -> argument; cat() logging dropped.
simulate_training_doublets <- function(counts, hc_idx, hc_types,
                                       types = sort(unique(hc_types)),
                                       n_per_pair = 200, seed = 1L) {
  set.seed(seed)
  type_to_cells <- split(seq_along(hc_idx), hc_types)
  het_pairs <- utils::combn(types, 2, simplify = FALSE)

  dblA <- integer(0); dblB <- integer(0); dbl_pl <- character(0)
  for (pr in het_pairs) {
    a <- pr[1]; b <- pr[2]
    if (length(type_to_cells[[a]]) < 1 || length(type_to_cells[[b]]) < 1) next
    ca <- hc_idx[sample(type_to_cells[[a]], n_per_pair, replace = TRUE)]
    cb <- hc_idx[sample(type_to_cells[[b]], n_per_pair, replace = TRUE)]
    dblA <- c(dblA, ca); dblB <- c(dblB, cb)
    dbl_pl <- c(dbl_pl, rep(pair_label(a, b), n_per_pair))
  }

  dbl_counts <- counts[, dblA, drop = FALSE] + counts[, dblB, drop = FALSE]
  colnames(dbl_counts) <- paste0("traindbl_", seq_len(ncol(dbl_counts)))

  list(counts = dbl_counts, pair_label = dbl_pl, idx_A = dblA, idx_B = dblB)
}
