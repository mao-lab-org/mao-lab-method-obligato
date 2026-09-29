# Internal construction of candidate pairs from up to five top-scoring types.

make_candidates <- function(top) {
  if (length(top) < 2L) {
    stop("At least two scored cell types are required for composition.", call. = FALSE)
  }
  rest <- seq.int(2L, length(top))
  flex_c1 <- rep(top[1L], length(rest))
  flex_c2 <- top[rest]
  if (length(top) >= 3L) {
    flex_c1 <- c(flex_c1, rep(top[2L], length(top) - 2L))
    flex_c2 <- c(flex_c2, top[seq.int(3L, length(top))])
  }
  idx <- utils::combn(length(top), 2L)
  list(
    fix1 = data.frame(c1 = rep(top[1L], length(rest)), c2 = top[rest],
                      stringsAsFactors = FALSE),
    flex2 = data.frame(c1 = flex_c1, c2 = flex_c2, stringsAsFactors = FALSE),
    all5 = data.frame(c1 = top[idx[1L, ]], c2 = top[idx[2L, ]],
                      stringsAsFactors = FALSE)
  )
}

# Internal composition ranker. For each cell, candidate pairs are formed from
# its five highest PhiSpace scores (or every available type when fewer than five
# exist) and ranked by their fitted logit-normal likelihood. `flagged` only
# selects which rows to return; it does not affect a cell's composition.
#
# `same_models` (optional) adds ONE same-type candidate per cell, the top-scoring
# type paired with itself ({top1, top1}); with the default all5 set that is 11
# candidates. `log_prior` (optional) is a named log prior over pair labels, added
# to each candidate's log-likelihood for RANKING only; the `ll` column of topk
# output stays the raw log-likelihood and a `log_prior` column is added. With
# both NULL the behaviour is exactly the heterotypic-only ranker.
compose_pairs <- function(S, pair_models, flagged = NULL,
                          candidates = c("all5", "flex2", "fix1"),
                          output = c("top1", "topk", "conformal"),
                          k = 3L, same_models = NULL, log_prior = NULL) {
  candidates <- match.arg(candidates)
  output <- match.arg(output)
  if (output == "conformal") {
    stop("Use compose(..., output = \"conformal\") for calibrated prediction sets.",
         call. = FALSE)
  }
  S <- .validate_scores(S, nrow(S), "S")
  if (!is.list(pair_models) || length(pair_models) < 1L ||
      is.null(names(pair_models)) || anyDuplicated(names(pair_models))) {
    stop("`pair_models` must be a non-empty, uniquely named list.", call. = FALSE)
  }
  if (!is.null(flagged)) {
    if (!is.logical(flagged) || length(flagged) != nrow(S) || anyNA(flagged)) {
      stop("`flagged` must be a non-missing logical vector with one value per cell.",
           call. = FALSE)
    }
  }
  k <- .validate_scalar_integer(k, "k", 1L)
  if (!is.null(same_models)) {
    if (!is.list(same_models) || !length(same_models) || is.null(names(same_models)) ||
        anyDuplicated(names(same_models))) {
      stop("`same_models` must be NULL or a non-empty, uniquely named list.", call. = FALSE)
    }
    sp <- strsplit(names(same_models), " + ", fixed = TRUE)
    if (!all(vapply(sp, function(p) length(p) == 2L && p[1L] == p[2L], logical(1)))) {
      stop("`same_models` must be named \"A + A\" (same-type pairs only).", call. = FALSE)
    }
    if (any(names(same_models) %in% names(pair_models))) {
      stop("`same_models` and `pair_models` must not share names.", call. = FALSE)
    }
  }
  all_models <- c(pair_models, same_models)
  if (!is.null(log_prior)) {
    if (!is.numeric(log_prior) || is.null(names(log_prior)) ||
        !all(names(all_models) %in% names(log_prior)) ||
        anyNA(log_prior[names(all_models)])) {
      stop("`log_prior` must be a named numeric vector covering every model.", call. = FALSE)
    }
  }

  n    <- nrow(S)
  dims <- colnames(S)
  Y    <- to_logit(S)
  idx  <- if (is.null(flagged)) seq_len(n) else which(flagged)
  n_top <- min(5L, ncol(S))
  top <- t(vapply(seq_len(n), function(i)
    dims[order(S[i, ], decreasing = TRUE)[seq_len(n_top)]],
    character(n_top)))

  pm_names <- names(all_models)

  # Candidate pairs are built from the SCORE COLUMN NAMES, so the pair models
  # must be keyed on the same cell-type vocabulary. If they are not, every
  # candidate misses and the function would return empty calls for every droplet
  # -- a silent wrong answer. Fail here instead, naming the cause.
  pm_types <- unique(unlist(strsplit(pm_names, " + ", fixed = TRUE)))
  if (!any(pm_types %in% dims)) {
    stop("No pair model matches the score columns: the models are named with ",
         "cell types (e.g. \"", pm_types[[1L]], "\") that do not appear among ",
         "the score columns (e.g. \"", dims[[1L]], "\"). The pair models and the ",
         "scores must use the same cell-type vocabulary; this usually means the ",
         "labels used to fit the models came from a different annotation than ",
         "the reference the scores were computed against.", call. = FALSE)
  }

  LL <- matrix(-Inf, n, length(pm_names), dimnames = list(NULL, pm_names))
  for (pl in pm_names) LL[, pl] <- mvn_ll_rows(Y, all_models[[pl]])

  rank_cell <- function(i) {
    cand <- make_candidates(top[i, ])[[candidates]]
    if (!is.null(same_models)) {
      # Appended last, so a tie with a heterotypic candidate goes to the latter.
      cand <- rbind(cand, data.frame(c1 = top[i, 1L], c2 = top[i, 1L],
                                     stringsAsFactors = FALSE))
    }
    pls  <- pair_label(cand$c1, cand$c2)
    lls  <- rep(-Inf, nrow(cand))
    in_m <- pls %in% pm_names
    if (any(in_m)) lls[in_m] <- LL[i, pls[in_m]]
    lpr  <- if (is.null(log_prior)) rep(0, length(pls)) else
      ifelse(in_m, unname(log_prior[pls]), 0)
    ord <- order(-(lls + lpr), seq_along(lls))
    list(cand = cand, ord = ord[in_m[ord]], lls = lls, lpr = lpr, has_model = any(in_m))
  }

  if (output == "top1") {
    c1 <- c2 <- rep(NA_character_, length(idx))
    for (m in seq_along(idx)) {
      r <- rank_cell(idx[m])
      if (r$has_model) {
        best <- r$ord[1L]
        c1[m] <- r$cand$c1[best]
        c2[m] <- r$cand$c2[best]
      }
    }
    return(data.frame(cell = idx, c1 = c1, c2 = c2, stringsAsFactors = FALSE))
  }

  empty <- data.frame(cell = integer(), rank = integer(), c1 = character(),
                      c2 = character(), ll = numeric(), stringsAsFactors = FALSE)
  if (!length(idx)) return(empty)
  out <- vector("list", length(idx))
  for (m in seq_along(idx)) {
    r <- rank_cell(idx[m])
    take <- if (r$has_model) r$ord[seq_len(min(k, length(r$ord)))] else integer()
    out[[m]] <- data.frame(
      cell = rep.int(idx[m], length(take)), rank = seq_along(take),
      c1 = r$cand$c1[take], c2 = r$cand$c2[take], ll = r$lls[take],
      stringsAsFactors = FALSE)
    if (!is.null(log_prior)) out[[m]]$log_prior <- r$lpr[take]
  }
  out <- Filter(NROW, out)
  if (!is.null(log_prior)) empty$log_prior <- numeric()
  if (!length(out)) empty else do.call(rbind, out)
}

#' Identify doublet composition with an optional detector
#'
#' Scores a count matrix against one or more labelled references, learns
#' logit-normal models for synthetic heterotypic and (by default) same-type
#' doublets, and assigns candidate cell-type pairs. The built-in detector combines PhiSpace score features with
#' library-size features. Composition is calculated for every query cell; use
#' `flag` to select calls made by the chosen detector.
#'
#' @param query A genes-by-cells matrix-like object containing raw non-negative counts.
#' @param reference A labelled `SingleCellExperiment`; for multiple annotation
#'   levels, a list with one reference per element of `phenotypes`.
#' @param phenotypes Character vector naming the label column in each reference.
#'   The first level defines composition; later levels add detection features.
#' @param hc_singlets Non-missing logical vector identifying high-confidence singlets.
#' @param hc_labels Optional labels for the high-confidence singlets. Supply either
#'   one value per query cell or one value per `TRUE` in `hc_singlets`.
#' @param detector `"builtin"`, `"none"`, a logical call vector, or a
#'   numeric score vector. Numeric scores require `detector_threshold`.
#' @param detector_threshold Threshold for a numeric external detector score.
#' @param n_per_pair Number of synthetic doublets generated per heterotypic pair.
#' @param candidates Candidate-pair strategy.
#' @param same_type Also consider a same-type composition (default `TRUE`): the
#'   top-scoring type paired with itself is added as one extra candidate per
#'   droplet, scored against a model fitted to simulated doublets of two cells of
#'   that type. `FALSE` restores the heterotypic-only composition exactly. Affects
#'   composition only; see "Same-type compositions and the pair prior".
#' @param pair_prior `"frequency"` (default) weights each candidate pair by how
#'   often two randomly drawn high-confidence singlets would form it,
#'   \eqn{\pi(\{a,b\}) = 2 p_a p_b} and \eqn{\pi(\{a,a\}) = p_a^2}; `"none"`
#'   ranks by likelihood alone. Affects composition only.
#' @param output One top pair, the top `k` pairs, or a conformal set.
#' @param k Number of pairs returned when `output = "topk"`.
#' @param alpha Miscoverage level for conformal sets.
#' @param seed Non-negative integer controlling simulation, scoring, and training.
#' @param eb Shrink the doublet-pair covariances toward a target predicted from
#'   the constituents' singlet models (default `TRUE`). `FALSE` restores the
#'   previous Ledoit-Wolf-only behaviour exactly.
#' @param eb_lambda A number in `[0, 1]` fixing the shrinkage weight (default
#'   `0.75`), or `"cv"` to choose it by cross-validated detection AUPRC on the
#'   simulated doublets. Cross-validation needs no ground truth and can return 0,
#'   declining to shrink; see the caution below on what it can resolve.
#' @param eb_grid Candidate weights for `eb_lambda = "cv"`; must contain 0.
#' @param eb_R Monte-Carlo draws used to build each target.
#' @param eb_nfold Cross-validation folds for `eb_lambda = "cv"`.
#' @param store_features Keep the per-droplet detection features in the result
#'   (default `TRUE`). The features are computed either way when the built-in
#'   detector is used; storing them costs roughly
#'   `(7 * levels + 3) * ncol(query) * 8` bytes, about 1.4 MB for 10,000
#'   droplets at two annotation levels and 140 MB at one million. Set `FALSE`
#'   for very large queries. With a supplied `detector` the features are not
#'   otherwise needed, so `TRUE` adds one pair-likelihood pass.
#' @param score_fn Advanced scoring function with the same interface as
#'   `score_phispace()`.
#' @return An `obligato_composition` object containing composition calls,
#'   detection scores and flags, the per-droplet detection features, fitted
#'   composition models, parameters, and provenance.
#'
#' @section Per-droplet detection features:
#'
#' `detection_features` holds one row per query droplet, row-named by droplet,
#' in query column order. Alongside the score summaries and the three
#' library-size features it carries `sing_LL`, `dbl_LL` and their difference
#' `ll_diff`: how much better the best two-type model explains a droplet than
#' the best one-type model. That quantity answers a question the library-size
#' statistic cannot, namely whether a droplet of entirely ordinary size is
#' doublet-like at all, and it is the basis for judging whether a group of calls
#' is supported.
#'
#' `ll_diff` compares the best pair over all pairs against the best singlet. To
#' assess the pair that was actually reported instead, subtract `sing_LL` from
#' the `ll` of that droplet's row in `composition`; the two differ whenever the
#' reported pair is not the global maximum, which the candidate rule does not
#' guarantee.
#'
#' With a supplied `detector`, only the primary annotation level contributes.
#' The per-level models exist solely inside the built-in detector's path, and
#' the composition models are keyed on the primary level's types, so a second
#' level's scores cannot be evaluated against them.
#'
#' @section Same-type compositions and the pair prior:
#'
#' A doublet can contain two cells of the same type. With `same_type = TRUE` each
#' droplet's candidates are the heterotypic pairs among its five highest-scoring
#' types plus one same-type candidate, the top type paired with itself. Its model
#' is fitted to simulated doublets built from two distinct high-confidence
#' singlets of that type. On experimentally captured doublets this names about
#' 55-64% of real same-type doublets correctly, at a cost of about 3 percentage
#' points on heterotypic doublets, almost all between closely related subtypes.
#'
#' The same-type training doublets are scored in their own pool with the
#' high-confidence singlets, separately from the heterotypic training pool, so
#' the detector, its features and its scores are identical whether or not
#' `same_type` is used.
#'
#' In score space a same-type doublet closely resembles a singlet of that type,
#' so a same-type call means "only this type is evident"; whether the droplet
#' holds one cell or two is the detector's question, not composition's.
#'
#' `pair_prior = "frequency"` adds the log prior to each candidate's
#' log-likelihood before ranking. When the scores cannot separate two
#' candidates, the pair made of more common types wins. This raises accuracy on
#' the realistic mix of doublets and lowers it when every pair is weighted
#' equally (balanced accuracy), because rarer-type pairs are pulled toward
#' common explanations. The frequencies come from the high-confidence singlets.
#' Conformal output (`output = "conformal"`) ignores both options.
#'
#' @section Shrinkage of the doublet-pair covariances:
#'
#' Each pair model needs a covariance over the cell-type score vector, and a
#' reference with `T` types has `T(T-1)/2` pairs to estimate, each from the
#' simulated doublets of that pair alone. When the score dimension is large
#' relative to those samples the estimates are noisy, and the noise is worst
#' exactly where the data are thinnest.
#'
#' With `eb = TRUE` (the default) each pair covariance is blended toward a
#' target predicted from the two constituents' *singlet* models,
#' \deqn{\Sigma_{eb} = (1-\lambda)\,\Sigma_{lw} + \lambda\,\Sigma_0,}
#' where \eqn{\Sigma_0} follows from the fact that a doublet's counts are the sum
#' of its constituents', so its score vector is approximately a depth-weighted
#' mixture of theirs. Because every type appears in many pairs, its singlet model
#' is estimated from far more cells than any single pair model, which is where
#' the strength being borrowed comes from. Singlet models are never shrunk.
#'
#' **The shrunk models are used for the detection features only.** Composition is
#' computed from the unshrunk models, so `compose()` returns identical
#' composition calls with and without `eb`. That split is empirical: across the
#' three benchmark datasets shrinkage improves detection at every weight tested
#' (+0.5 to +1.3 percentage points of AUPRC) and does not improve deconvolution
#' (-0.31, -0.03 and 0.00 points, none significant). The mechanism is that
#' shrinking a pair toward a target built from its own constituents stabilises
#' the *level* of the pair likelihood, which detection consumes as scalar
#' features, while making different pairs more alike — which is precisely what
#' the arg-max across pairs has to separate.
#'
#' The default weight is fixed at `0.75` rather than cross-validated. The
#' detection response is flat between about 0.7 and 0.9, so a fixed value costs
#' little against per-dataset tuning, and 0.75 keeps a margin from
#' \eqn{\lambda = 1}, where the pair's own simulated doublets stop contributing
#' altogether and within-lineage detection degrades sharply (up to -8 points of
#' AUPRC on one benchmark dataset).
#'
#' `eb_lambda = "cv"` selects the weight by cross-validated **detection AUPRC**,
#' matching the task the shrinkage is used for. Be aware of its resolution: it
#' detects reliably that shrinkage helps — the deficit at \eqn{\lambda = 0} is
#' around four times the fold-to-fold standard deviation — but across
#' \eqn{\lambda \ge 0.5} the spread is about the size of that standard
#' deviation, so the selected value in that region is close to arbitrary and can
#' land on the endpoint. It is offered as an option, not as the default, and it
#' costs a detector refit per fold per grid point.
#'
#' `eb_grid` must contain 0, since that is what lets the procedure decline to
#' shrink. `eb = FALSE` restores the previous Ledoit-Wolf-only behaviour exactly.
#'
#' The chosen weight and the cross-validation curve are recorded in `info$eb`.
#'
#' @export
compose <- function(query, reference, phenotypes, hc_singlets,
                    hc_labels = NULL, detector = "builtin",
                    detector_threshold = NULL, n_per_pair = 200L,
                    candidates = c("all5", "flex2", "fix1"),
                    same_type = TRUE, pair_prior = c("frequency", "none"),
                    output = c("top1", "topk", "conformal"),
                    k = 3L, alpha = 0.10, seed = 1L,
                    eb = TRUE, eb_lambda = 0.75,
                    eb_grid = seq(0, 1, by = 0.1), eb_R = 2000L,
                    eb_nfold = 5L, store_features = TRUE,
                    score_fn = score_phispace) {
  .validate_counts(query)
  candidates <- match.arg(candidates)
  pair_prior <- match.arg(pair_prior)
  output <- match.arg(output)
  if (!is.logical(same_type) || length(same_type) != 1L || is.na(same_type)) {
    stop("`same_type` must be TRUE or FALSE.", call. = FALSE)
  }
  n_per_pair <- .validate_scalar_integer(n_per_pair, "n_per_pair", 30L)
  seed <- .validate_scalar_integer(seed, "seed", 0L)
  k <- .validate_scalar_integer(k, "k", 1L)
  if (!is.logical(store_features) || length(store_features) != 1L ||
      is.na(store_features)) {
    stop("`store_features` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(hc_singlets) || length(hc_singlets) != ncol(query) ||
      anyNA(hc_singlets)) {
    stop("`hc_singlets` must be a non-missing logical vector with one value per cell.",
         call. = FALSE)
  }
  if (!length(hc_singlets) || !any(hc_singlets)) {
    stop("`hc_singlets` must select at least one cell.", call. = FALSE)
  }
  if (output == "conformal" &&
      (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) ||
       alpha <= 0 || alpha >= 1)) {
    stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  cnt <- query
  hc_idx <- which(hc_singlets)
  phenotypes <- as.character(phenotypes)
  if (!length(phenotypes) || anyNA(phenotypes) || any(!nzchar(phenotypes))) {
    stop("`phenotypes` must contain one or more non-empty column names.", call. = FALSE)
  }
  refs <- if (length(phenotypes) == 1L) list(reference) else reference
  if (!is.list(refs) || length(refs) != length(phenotypes)) {
    stop("For multiple `phenotypes`, `reference` must be a list of equal length.",
         call. = FALSE)
  }
  labels_by_level <- lapply(seq_along(refs), function(i) {
    z <- tryCatch(refs[[i]][[phenotypes[[i]]]], error = function(e) NULL)
    if (is.null(z) || !length(z)) {
      stop(sprintf("Reference %d has no `%s` annotation.", i, phenotypes[[i]]),
           call. = FALSE)
    }
    as.character(z)
  })
  ref1 <- refs[[1L]]
  pheno1 <- phenotypes[[1L]]

  S_all <- .with_seed(seed,
    score_fn(cnt, ref1, pheno1, " all-cells"))
  S_all <- .validate_scores(S_all, ncol(cnt), "Primary query scores")

  if (is.null(hc_labels)) {
    hc_labels <- colnames(S_all)[max.col(
      S_all[hc_idx, , drop = FALSE], ties.method = "first")]
  } else if (length(hc_labels) == ncol(cnt)) {
    hc_labels <- hc_labels[hc_idx]
  } else if (length(hc_labels) != length(hc_idx)) {
    stop("`hc_labels` must have one value per query cell or selected singlet.",
         call. = FALSE)
  }
  hc_labels <- as.character(hc_labels)
  ref_types <- sort(unique(labels_by_level[[1L]]))
  valid <- !is.na(hc_labels) & hc_labels %in% ref_types
  if (any(!valid)) {
    warning(sprintf("Dropping %d high-confidence singlet(s) with missing or unknown labels.",
                    sum(!valid)), call. = FALSE)
  }
  hc_idx <- hc_idx[valid]
  hc_labels <- hc_labels[valid]
  types <- sort(unique(hc_labels))
  if (length(types) < 2L) {
    stop("At least two valid high-confidence singlet types are required.", call. = FALSE)
  }

  sim <- simulate_training_doublets(cnt, hc_idx, hc_labels, types,
                                    n_per_pair, seed)
  n_sing <- length(hc_idx)
  joint_counts <- cbind(cnt[, hc_idx, drop = FALSE], sim$counts)
  S_join <- .with_seed(seed + 1L,
    score_fn(joint_counts, ref1, pheno1, " train-joint"))
  S_join <- .validate_scores(S_join, ncol(joint_counts), "Primary training scores",
                             expected = colnames(S_all))
  S_sing <- S_join[seq_len(n_sing), , drop = FALSE]
  S_dbl <- S_join[n_sing + seq_len(ncol(sim$counts)), , drop = FALSE]

  mods <- fit_doublet_models(S_sing, hc_labels, S_dbl, sim$pair_label, types)

  # Same-type composition models. Simulated from two distinct singlets of one
  # type and scored in their OWN pool with the singlets, so the heterotypic pool
  # above -- and with it every detection feature and score -- is untouched.
  same_models <- NULL
  if (isTRUE(same_type)) {
    sim_same <- simulate_same_type_doublets(cnt, hc_idx, hc_labels, types,
                                            n_per_pair, seed + 2L)
    if (ncol(sim_same$counts)) {
      same_counts <- cbind(cnt[, hc_idx, drop = FALSE], sim_same$counts)
      S_same_join <- .with_seed(seed + 3L,
        score_fn(same_counts, ref1, pheno1, " train-same"))
      S_same_join <- .validate_scores(S_same_join, ncol(same_counts),
                                      "Same-type training scores",
                                      expected = colnames(S_all))
      same_models <- fit_doublet_models(
        S_same_join[seq_len(n_sing), , drop = FALSE], hc_labels,
        S_same_join[n_sing + seq_len(ncol(sim_same$counts)), , drop = FALSE],
        sim_same$pair_label, types)$pair_models
      if (!length(same_models)) same_models <- NULL
    }
  }
  type_freq <- as.numeric(table(factor(hc_labels, levels = types))) / length(hc_labels)
  names(type_freq) <- types
  log_prior <- if (pair_prior == "frequency")
    pair_log_prior(c(names(mods$pair_models), names(same_models)), type_freq) else NULL

  # Empirical-Bayes shrinkage of the PAIR covariances is applied to the DETECTION
  # path only, and is set up inside the detector branch below. `mods` stays
  # UNSHRUNK here and is what compose_pairs() and the returned object use, so
  # composition output is identical with and without `eb`.
  #
  # Measured on the three benchmark datasets: shrinkage improves detection at
  # every step of lambda (+0.5 to +1.3 pp AUPRC) and does not improve
  # deconvolution (-0.31 / -0.03 / 0.00 pp, none significant). Shrinking toward a
  # target built from a pair's own constituents stabilises the LEVEL of the pair
  # likelihood -- which detection consumes as scalar features -- while making
  # pairs more alike, which is what the arg-max across pairs must separate.
  eb_info <- NULL
  if (!identical(eb_lambda, "cv") &&
      (!is.numeric(eb_lambda) || length(eb_lambda) != 1L || is.na(eb_lambda) ||
       eb_lambda < 0 || eb_lambda > 1)) {
    stop("`eb_lambda` must be \"cv\" or a single number in [0, 1].", call. = FALSE)
  }
  if (!0 %in% eb_grid) {
    stop("`eb_grid` must contain 0, so that shrinkage can be declined.", call. = FALSE)
  }
  # Only warn if the weight was actually set: the default is now a number, so
  # comparing against "cv" would fire for anyone who merely passed eb = FALSE.
  if (!isTRUE(eb) && !identical(eb_lambda, .EB_LAMBDA_DEFAULT)) {
    warning("`eb_lambda` is ignored when `eb = FALSE`.", call. = FALSE)
  }
  if (!length(mods$pair_models)) {
    stop("No doublet-pair models could be fitted; increase `n_per_pair`.", call. = FALSE)
  }

  det_score <- threshold <- flag <- detection_model <- NULL
  # Per-droplet detection features. `ll_diff` here is the difference between the
  # best pair log-likelihood and the best singlet log-likelihood, i.e. how much
  # better a two-type model explains the droplet than any one-type model. It is
  # computed for the detector and was previously discarded; retaining it is what
  # lets a caller ask whether a group of calls is doublet-like at all, which the
  # library-size statistic cannot answer for a droplet of ordinary size.
  F_all <- NULL

  if (identical(detector, "builtin")) {
    if (!length(mods$sing_models)) {
      stop("No singlet models could be fitted; provide at least 15 cells for one type.",
           call. = FALSE)
    }
    n_levels <- length(phenotypes)
    suffix <- function(i) if (n_levels > 1L) paste0("_l", i) else ""

    # Score every annotation level ONCE. Model fitting is separated from scoring
    # so that the cross-validation below can refit models per fold without
    # re-scoring, which is the expensive part.
    level_scores <- function(i) {
      if (i == 1L) return(list(Sa = S_all, Ss = S_sing, Sd = S_dbl))
      Sa <- .with_seed(seed + 100L + i,
        score_fn(cnt, refs[[i]], phenotypes[[i]], sprintf(" all-cells L%d", i)))
      Sa <- .validate_scores(Sa, ncol(cnt), sprintf("Query scores for level %d", i))
      Sj <- .with_seed(seed + 200L + i,
        score_fn(joint_counts, refs[[i]], phenotypes[[i]],
                 sprintf(" train-joint L%d", i)))
      Sj <- .validate_scores(Sj, ncol(joint_counts),
                             sprintf("Training scores for level %d", i),
                             expected = colnames(Sa))
      list(Sa = Sa,
           Ss = Sj[seq_len(n_sing), , drop = FALSE],
           Sd = Sj[n_sing + seq_len(ncol(sim$counts)), , drop = FALSE])
    }
    lvl <- lapply(seq_len(n_levels), level_scores)

    # Select the shrinkage weight for the detection features.
    eb_lam <- 0
    if (isTRUE(eb)) {
      lib_sing <- as.numeric(Matrix::colSums(cnt[, hc_idx, drop = FALSE]))
      if (identical(eb_lambda, "cv")) {
        sel <- .cv_lambda_detection(lvl, hc_labels, sim$pair_label, types,
                                    cnt[, hc_idx, drop = FALSE], sim$counts,
                                    lib_sing, grid = eb_grid, R = eb_R,
                                    nfold = eb_nfold, seed = seed,
                                    suffix = suffix)
        eb_lam <- sel$lambda
        eb_info <- list(lambda = eb_lam, cv = sel$cv, grid = eb_grid,
                        requested = "cv", criterion = "detection_auprc",
                        applies_to = "detection", n_shrunk = NA_integer_)
      } else {
        eb_lam <- eb_lambda
        eb_info <- list(lambda = eb_lam, cv = NULL, grid = eb_grid,
                        requested = eb_lambda, criterion = "fixed",
                        applies_to = "detection", n_shrunk = NA_integer_)
      }
    }

    level_features <- function(i) {
      L <- lvl[[i]]
      m <- if (i == 1L) mods else
        fit_doublet_models(L$Ss, hc_labels, L$Sd, sim$pair_label, types)
      if (eb_lam > 0) {
        fit <- apply_eb_pairs(m, L$Ss, hc_labels, L$Sd, sim$pair_label,
                              as.numeric(Matrix::colSums(cnt[, hc_idx, drop = FALSE])),
                              lambda = eb_lam, grid = eb_grid,
                              R = eb_R, seed = seed)
        m <- fit$models
        if (i == 1L) eb_info$n_shrunk <<- fit$n_shrunk
      }
      if (!length(m$sing_models) || !length(m$pair_models)) {
        stop(sprintf("No usable models could be fitted for annotation level %d.", i),
             call. = FALSE)
      }
      list(
        sing = compute_features(L$Ss, m$sing_models, m$pair_models, suffix(i)),
        dbl = compute_features(L$Sd, m$sing_models, m$pair_models, suffix(i)),
        all = compute_features(L$Sa, m$sing_models, m$pair_models, suffix(i)),
        models = m)
    }
    parts <- lapply(seq_len(n_levels), level_features)
    F_sing <- do.call(cbind, lapply(parts, `[[`, "sing"))
    F_dbl <- do.call(cbind, lapply(parts, `[[`, "dbl"))
    F_all <- do.call(cbind, lapply(parts, `[[`, "all"))

    cnt_sing <- cnt[, hc_idx, drop = FALSE]
    type_median <- libsize_type_medians(cnt_sing, S_sing)
    F_sing <- cbind(F_sing, libsize_features(
      cnt_sing, top1_type(S_sing), type_median))
    F_dbl <- cbind(F_dbl, libsize_features(
      sim$counts, top1_type(S_dbl), type_median))
    F_all <- cbind(F_all, libsize_features(
      cnt, top1_type(S_all), type_median))

    bst <- train_detector(F_sing, F_dbl, seed = seed)
    det_score <- as.numeric(predict(bst, as.matrix(F_all)))
    threshold <- ecdf_threshold(
      det_score, as.numeric(predict(bst, as.matrix(F_dbl))))
    flag <- det_score > threshold
    detection_model <- list(
      classifier = bst, threshold = threshold, feature_names = names(F_all),
      type_median = type_median, phenotypes = phenotypes,
      level_models = lapply(parts, `[[`, "models"))
  } else if (is.logical(detector)) {
    if (length(detector) != ncol(cnt) || anyNA(detector)) {
      stop("A logical `detector` must contain one non-missing call per cell.",
           call. = FALSE)
    }
    flag <- detector
  } else if (is.numeric(detector)) {
    if (length(detector) != ncol(cnt) || anyNA(detector) ||
        any(!is.finite(detector))) {
      stop("A numeric `detector` must contain one finite score per cell.",
           call. = FALSE)
    }
    if (length(detector_threshold) != 1L || !is.numeric(detector_threshold) ||
        is.na(detector_threshold) || !is.finite(detector_threshold)) {
      stop("Numeric detector scores require a finite `detector_threshold`.",
           call. = FALSE)
    }
    det_score <- as.numeric(detector)
    threshold <- detector_threshold
    flag <- det_score > threshold
  } else if (!identical(detector, "none")) {
    stop("`detector` must be \"builtin\", \"none\", or a logical/numeric vector.",
         call. = FALSE)
  }

  # Supplied detector: the built-in branch above did not run, so the features do
  # not exist yet. `mods` is fitted regardless because composition needs it, so
  # they can still be produced -- at one extra pair-likelihood pass, which is why
  # it is gated on the caller asking. Only the PRIMARY annotation level is
  # available here: the per-level models and score matrices live inside the
  # built-in branch, and `mods` is keyed on the primary level's types, so a
  # second level's scores would be evaluated against the wrong label space.
  if (is.null(F_all) && isTRUE(store_features) && !identical(detector, "none")) {
    F_all <- cbind(
      compute_features(S_all, mods$sing_models, mods$pair_models, ""),
      libsize_features(cnt, top1_type(S_all),
                       libsize_type_medians(cnt[, hc_idx, drop = FALSE], S_sing)))
  }

  conformal <- NULL
  if (output == "conformal") {
    sim_cal <- simulate_training_doublets(
      cnt, hc_idx, hc_labels, types, n_per_pair, seed + 1L)
    cal_counts <- cbind(cnt, sim_cal$counts)
    S_cal_join <- .with_seed(seed + 1000L,
      score_fn(cal_counts, ref1, pheno1, " conformal-cal"))
    S_cal_join <- .validate_scores(
      S_cal_join, ncol(cal_counts), "Conformal scores",
      expected = colnames(S_all))
    S_cal_query <- S_cal_join[seq_len(ncol(cnt)), , drop = FALSE]
    S_cal <- S_cal_join[ncol(cnt) + seq_len(ncol(sim_cal$counts)), , drop = FALSE]
    log_h_cal <- all_pair_ll(S_cal, mods$pair_models)
    log_h_query <- all_pair_ll(S_cal_query, mods$pair_models)
    cf <- conformal_pairs(log_h_cal, sim_cal$pair_label, log_h_query, alpha)

    rows <- lapply(seq_len(ncol(cnt)), function(i) {
      ps <- cf$sets[[i]]
      if (!length(ps)) return(NULL)
      pair_parts <- strsplit(ps, " \\+ ")
      data.frame(
        cell = rep.int(i, length(ps)),
        c1 = vapply(pair_parts, `[[`, character(1), 1L),
        c2 = vapply(pair_parts, `[[`, character(1), 2L),
        pair = ps, ll = unname(log_h_query[i, ps]),
        stringsAsFactors = FALSE)
    })
    rows <- Filter(NROW, rows)
    comp <- if (length(rows)) do.call(rbind, rows) else
      data.frame(cell = integer(), c1 = character(), c2 = character(),
                 pair = character(), ll = numeric(), stringsAsFactors = FALSE)
    conformal <- list(q_hat = cf$q_hat, alpha = alpha, sizes = cf$sizes)
  } else {
    comp <- compose_pairs(S_all, mods$pair_models, candidates = candidates,
                          output = output, k = k, same_models = same_models,
                          log_prior = log_prior)
  }

  # One row per query droplet, in query column order. `flag` and
  # `detection_score` are unnamed vectors, so naming these rows is what makes a
  # join to the composition table unambiguous rather than positional.
  if (!is.null(F_all)) {
    if (nrow(F_all) != ncol(cnt)) {
      stop("Internal error: detection features have ", nrow(F_all),
           " rows for ", ncol(cnt), " droplets.", call. = FALSE)
    }
    rownames(F_all) <- colnames(cnt)
  }
  if (!isTRUE(store_features)) F_all <- NULL

  detector_name <- if (is.character(detector)) detector else "byo"
  info <- list(
    eb = eb_info,
    n_cells = ncol(cnt), n_hc_singlets = n_sing, types = types,
    phenotypes = phenotypes, detector = detector_name,
    detection_feature_source = if (identical(detector, "builtin"))
      "PhiSpace scores plus library size" else NULL,
    candidates = candidates, same_type = same_type,
    n_same_models = length(same_models), pair_prior = pair_prior,
    type_freq = type_freq, output = output, n_per_pair = n_per_pair,
    seed = seed, threshold = threshold,
    alpha = if (output == "conformal") alpha else NULL)

  structure(
    list(composition = comp, detection_score = det_score, flag = flag,
         detection_features = F_all,
         conformal = conformal, models = c(mods, list(same_models = same_models)),
         detection_model = detection_model, info = info,
         provenance = .package_provenance()),
    class = "obligato_composition")
}
