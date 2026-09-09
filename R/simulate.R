# Internal synthetic training-doublet generator. For each heterotypic type pair,
# it samples high-confidence singlets with replacement and sums their raw counts.
# Sampling is reproducible without modifying the caller's RNG state.
simulate_training_doublets <- function(counts, hc_idx, hc_types,
                                       types = sort(unique(hc_types)),
                                       n_per_pair = 200L, seed = 1L) {
  if (length(hc_idx) != length(hc_types) || anyNA(hc_idx) || anyNA(hc_types)) {
    stop("`hc_idx` and `hc_types` must be complete vectors of equal length.",
         call. = FALSE)
  }
  if (length(types) < 2L) {
    stop("At least two singlet types are required to simulate doublets.",
         call. = FALSE)
  }
  n_per_pair <- .validate_scalar_integer(n_per_pair, "n_per_pair", 1L)
  .with_seed(seed, {
    type_to_cells <- split(seq_along(hc_idx), hc_types)
    het_pairs <- utils::combn(types, 2L, simplify = FALSE)

    dblA <- integer(); dblB <- integer(); dbl_pl <- character()
    for (pr in het_pairs) {
      a <- pr[1L]; b <- pr[2L]
      if (!length(type_to_cells[[a]]) || !length(type_to_cells[[b]])) next
      ca <- hc_idx[sample(type_to_cells[[a]], n_per_pair, replace = TRUE)]
      cb <- hc_idx[sample(type_to_cells[[b]], n_per_pair, replace = TRUE)]
      dblA <- c(dblA, ca); dblB <- c(dblB, cb)
      dbl_pl <- c(dbl_pl, rep(pair_label(a, b), n_per_pair))
    }

    dbl_counts <- counts[, dblA, drop = FALSE] + counts[, dblB, drop = FALSE]
    colnames(dbl_counts) <- paste0("traindbl_", seq_len(ncol(dbl_counts)))
    list(counts = dbl_counts, pair_label = dbl_pl, idx_A = dblA, idx_B = dblB)
  })
}
