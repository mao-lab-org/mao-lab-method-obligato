##############################################################################
# Composition candidate construction.
# Extracted VERBATIM from run_methodC_compose.R (M1).
##############################################################################

make_candidates <- function(top5) {
  list(
    fix1 = data.frame(c1 = rep(top5[1], 4), c2 = top5[2:5], stringsAsFactors = FALSE),
    flex2 = data.frame(
      c1 = c(rep(top5[1], 4), rep(top5[2], 3)),
      c2 = c(top5[2:5], top5[c(3, 4, 5)]),
      stringsAsFactors = FALSE),
    all5 = {
      idx <- combn(5, 2)
      data.frame(c1 = top5[idx[1, ]], c2 = top5[idx[2, ]], stringsAsFactors = FALSE)
    }
  )
}

##############################################################################
# DblLik composition: which two cell types a doublet contains (M1b).
#
# Extracted from run_compose_kingetal.R:342-397 (Section 6), the DblLik branch,
# identical to run_methodC_compose.R's process_batch composition block.
#
# For each cell: take the top-5 scoring types, form candidate pairs for the chosen
# candidate set, score each pair under its fitted per-pair logit-normal model, and
# rank by likelihood. `output = "top1"` returns the single max-likelihood pair —
# the stored `dbllik_<candidates>_c1/c2` and the value M2b reproduces; `"topk"`
# returns the k best pairs per cell. `"conformal"` is NOT handled here (this is the
# low-level per-cell ranker with no calibration set); compose() intercepts
# output = "conformal" and runs split-conformal — see R/conformal.R. Calling
# compose_pairs directly with "conformal" errors.
#
# BRING YOUR OWN DETECTOR: `flagged` is any logical vector over the rows of `S`.
# When NULL every cell is composed (matching the script, which composes all cells
# then filters downstream); when supplied only the flagged cells are composed. The
# per-cell composition is identical either way — flagging only selects rows.
#
# Faithful-but-parameterised: `candidates` defaults to the locked "all5". The
# per-cell strict-argmax ("first candidate wins ties, models absent from the fit
# score -Inf") is preserved exactly via order(-ll, seq_along(ll)); the per-pair
# likelihoods are computed with the same mvn_ll_rows, batched over cells for speed
# (row-wise identical to the script's one-row-at-a-time calls).
compose_pairs <- function(S, pair_models, flagged = NULL,
                    candidates = c("all5", "flex2", "fix1"),
                    output = c("top1", "topk", "conformal"),
                    k = 3) {
  candidates <- match.arg(candidates)
  output <- match.arg(output)
  if (output == "conformal") {
    stop("output = \"conformal\" is not yet implemented; conformal calibration is a ",
         "separate step (package_plan.md Phase B). Use \"top1\" or \"topk\".",
         call. = FALSE)
  }

  n    <- nrow(S)
  dims <- colnames(S)
  Y    <- to_logit(S)
  idx  <- if (is.null(flagged)) seq_len(n) else which(flagged)

  # Top-5 scoring types per cell (identical order() call to the script).
  top5 <- t(apply(S, 1, function(r) dims[order(r, decreasing = TRUE)[1:5]]))

  # Likelihood of every cell under every fitted pair model, batched.
  pm_names <- names(pair_models)
  LL <- matrix(-Inf, n, length(pm_names), dimnames = list(NULL, pm_names))
  for (pl in pm_names) LL[, pl] <- mvn_ll_rows(Y, pair_models[[pl]])

  rank_cell <- function(i) {
    cand <- make_candidates(top5[i, ])[[candidates]]
    pls  <- pair_label(cand$c1, cand$c2)
    lls  <- rep(-Inf, nrow(cand))
    in_m <- pls %in% pm_names
    if (any(in_m)) lls[in_m] <- LL[i, pls[in_m]]
    ord  <- order(-lls, seq_along(lls))   # ties -> first candidate, matches script
    list(cand = cand, ord = ord, lls = lls)
  }

  if (output == "top1") {
    c1 <- character(length(idx)); c2 <- character(length(idx))
    for (m in seq_along(idx)) {
      r <- rank_cell(idx[m]); best <- r$ord[1]
      c1[m] <- r$cand$c1[best]; c2[m] <- r$cand$c2[best]
    }
    return(data.frame(cell = idx, c1 = c1, c2 = c2, stringsAsFactors = FALSE))
  }

  # output == "topk"
  out <- vector("list", length(idx))
  for (m in seq_along(idx)) {
    r <- rank_cell(idx[m])
    take <- r$ord[seq_len(min(k, length(r$ord)))]
    out[[m]] <- data.frame(
      cell = idx[m], rank = seq_along(take),
      c1 = r$cand$c1[take], c2 = r$cand$c2[take], ll = r$lls[take],
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

##############################################################################
# compose(): the package's main entry point (M1b orchestrator).
#
# Composition is the obligato part; detection is an OPTIONAL module. Given a query
# and a labelled reference, compose() scores the cells, simulates training
# doublets, fits the per-pair models, optionally runs detection, and returns the
# per-cell composition. It is a thin wrapper over the tested building blocks
# (simulate_training_doublets, score_phispace, fit_doublet_models, compute_features,
# train_detector, ecdf_threshold, compose_pairs); no method logic lives here.
#
# `detector` selects the detection module:
#   "builtin" (default) : train and apply our xgboost detector; return score + flag
#   <logical/score vec> : BRING YOUR OWN DETECTOR — use these flags, skip ours
#   "none"              : compose only, no flagging
#
# `hc_labels` may be supplied (e.g. the King precomputed PhiSpace top-1) or derived
# from the query scores. `score_fn` is injectable so the orchestration is testable
# without PhiSpace; it defaults to score_phispace and must return a cells x types
# score matrix given (counts, reference, phenotypes, label).
#
# REPRODUCIBILITY: deterministic given `seed` (see simulate_training_doublets and
# train_detector). It reproduces itself exactly; it does NOT bit-reproduce the old
# non-deterministic stored numbers — that comparison is tolerance-based (M2b).
compose <- function(query, reference, phenotypes, hc_singlets,
                    hc_labels = NULL,
                    detector = "builtin",
                    n_per_pair = 200,
                    candidates = c("all5", "flex2", "fix1"),
                    output = c("top1", "topk", "conformal"),
                    k = 3, alpha = 0.10, seed = 1L,
                    score_fn = score_phispace) {
  candidates <- match.arg(candidates)
  output     <- match.arg(output)
  stopifnot(is.logical(hc_singlets), length(hc_singlets) == ncol(query))

  cnt    <- query
  hc_idx <- which(hc_singlets)

  # Score all query cells (inference scores, and the source of derived hc labels).
  S_all <- score_fn(cnt, reference, phenotypes, " all-cells")

  # High-confidence singlet labels: supplied, or top-1 of their scores.
  if (is.null(hc_labels)) {
    hc_labels <- colnames(S_all)[max.col(S_all[hc_idx, , drop = FALSE], ties.method = "first")]
  }
  ref_types <- sort(unique(as.character(reference[[phenotypes]])))
  valid     <- hc_labels %in% ref_types
  hc_idx    <- hc_idx[valid]
  hc_labels <- hc_labels[valid]
  types     <- sort(unique(hc_labels))

  # Simulate training doublets, then score the joint pool (locked joint norm).
  sim   <- simulate_training_doublets(cnt, hc_idx, hc_labels, types, n_per_pair, seed)
  n_sing <- length(hc_idx)
  S_join <- score_fn(cbind(cnt[, hc_idx, drop = FALSE], sim$counts),
                     reference, phenotypes, " train-joint")
  S_sing <- S_join[seq_len(n_sing), , drop = FALSE]
  S_dbl  <- S_join[(n_sing + 1L):nrow(S_join), , drop = FALSE]

  mods <- fit_doublet_models(S_sing, hc_labels, S_dbl, sim$pair_label, types)

  # Optional detection module.
  det_score <- NULL; threshold <- NULL; flag <- NULL
  if (identical(detector, "builtin")) {
    F_sing <- compute_features(S_sing, mods$sing_models, mods$pair_models)
    F_dbl  <- compute_features(S_dbl,  mods$sing_models, mods$pair_models)
    F_all  <- compute_features(S_all,  mods$sing_models, mods$pair_models)
    bst       <- train_detector(F_sing, F_dbl, seed = seed)
    det_score <- predict(bst, as.matrix(F_all))
    threshold <- ecdf_threshold(det_score, predict(bst, as.matrix(F_dbl)))
    flag      <- det_score > threshold
  } else if (is.logical(detector) || is.numeric(detector)) {
    stopifnot(length(detector) == ncol(query))
    flag <- if (is.numeric(detector)) detector > 0.5 else detector
  } else if (!identical(detector, "none")) {
    stop("`detector` must be \"builtin\", \"none\", or a logical/numeric vector.", call. = FALSE)
  }

  conformal <- NULL
  if (output == "conformal") {
    # Split-conformal prediction sets. A real query has no ground-truth pairs, so
    # calibrate on a FRESH, independent batch of simulated doublets (known pairs,
    # drawn from the same per-pair models, distinct seed so it is not the fit set).
    # The calibration doublets are scored jointly with the query cells so they share
    # the query normalisation regime, then scored against the fitted pair models.
    sim_cal <- simulate_training_doublets(cnt, hc_idx, hc_labels, types,
                                          n_per_pair, seed + 1L)
    S_cal_join <- score_fn(cbind(cnt, sim_cal$counts), reference, phenotypes,
                           " conformal-cal")
    S_cal       <- S_cal_join[(ncol(cnt) + 1L):nrow(S_cal_join), , drop = FALSE]
    log_h_cal   <- all_pair_ll(S_cal, mods$pair_models)
    log_h_query <- all_pair_ll(S_all, mods$pair_models)

    cf <- conformal_pairs(log_h_cal, sim_cal$pair_label, log_h_query, alpha = alpha)

    # Long-format composition: one row per (cell, pair) in the cell's set. Cells
    # with an empty set contribute no rows (recover them from `conformal$sizes`).
    rows <- lapply(seq_len(nrow(S_all)), function(i) {
      ps <- cf$sets[[i]]
      if (length(ps) == 0L) return(NULL)
      parts <- strsplit(ps, " \\+ ")
      data.frame(cell = i,
                 c1  = vapply(parts, function(p) p[1], character(1)),
                 c2  = vapply(parts, function(p) p[2], character(1)),
                 pair = ps, ll = log_h_query[i, ps],
                 stringsAsFactors = FALSE)
    })
    comp      <- do.call(rbind, rows)
    conformal <- list(q_hat = cf$q_hat, alpha = alpha, sizes = cf$sizes)
  } else {
    comp <- compose_pairs(S_all, mods$pair_models, candidates = candidates,
                          output = output, k = k)
  }

  info <- list(n_cells = ncol(query), n_hc_singlets = n_sing, types = types,
               detector = if (is.character(detector)) detector else "byo",
               candidates = candidates, output = output,
               n_per_pair = n_per_pair, seed = seed, threshold = threshold,
               alpha = if (output == "conformal") alpha else NULL)

  structure(list(composition = comp, detection_score = det_score, flag = flag,
                 conformal = conformal, models = mods, info = info),
            class = "obligato_composition")
}
