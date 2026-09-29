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

# Internal synthetic SAME-TYPE doublet generator: for each type, sums the raw
# counts of two DISTINCT high-confidence singlets of that type. Used for the
# composition models only; the detector never sees these doublets. Types with
# fewer than two singlets are skipped. Reproducible without touching the caller's
# RNG state, and drawn from its own seed so the heterotypic stream is unchanged.
simulate_same_type_doublets <- function(counts, hc_idx, hc_types,
                                        types = sort(unique(hc_types)),
                                        n_per_pair = 200L, seed = 1L) {
  if (length(hc_idx) != length(hc_types) || anyNA(hc_idx) || anyNA(hc_types)) {
    stop("`hc_idx` and `hc_types` must be complete vectors of equal length.",
         call. = FALSE)
  }
  n_per_pair <- .validate_scalar_integer(n_per_pair, "n_per_pair", 1L)
  .with_seed(seed, {
    type_to_cells <- split(seq_along(hc_idx), hc_types)
    dblA <- integer(); dblB <- integer(); dbl_pl <- character()
    for (a in types) {
      pool <- type_to_cells[[a]]
      if (length(pool) < 2L) next
      ia <- sample(pool, n_per_pair, replace = TRUE)
      ib <- sample(pool, n_per_pair, replace = TRUE)
      while (any(clash <- ia == ib)) ib[clash] <- sample(pool, sum(clash), replace = TRUE)
      dblA <- c(dblA, hc_idx[ia]); dblB <- c(dblB, hc_idx[ib])
      dbl_pl <- c(dbl_pl, rep(pair_label(a, a), n_per_pair))
    }
    dbl_counts <- counts[, dblA, drop = FALSE] + counts[, dblB, drop = FALSE]
    colnames(dbl_counts) <- paste0("traindbl_same_", seq_len(ncol(dbl_counts)))
    list(counts = dbl_counts, pair_label = dbl_pl, idx_A = dblA, idx_B = dblB)
  })
}

# Internal log prior over unordered type pairs from type frequencies p:
# pi({a,b}) = 2 p_a p_b for a != b, pi({a,a}) = p_a^2 -- the probability that two
# independently drawn cells form that pair. Returned for the requested labels.
pair_log_prior <- function(labels, type_freq) {
  parts <- strsplit(labels, " + ", fixed = TRUE)
  out <- vapply(parts, function(p) {
    if (length(p) != 2L || !all(p %in% names(type_freq))) return(NA_real_)
    if (p[1L] == p[2L]) 2 * log(type_freq[[p[1L]]]) else
      log(2) + log(type_freq[[p[1L]]]) + log(type_freq[[p[2L]]])
  }, numeric(1))
  names(out) <- labels
  out
}
