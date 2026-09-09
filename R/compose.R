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
compose_pairs <- function(S, pair_models, flagged = NULL,
                          candidates = c("all5", "flex2", "fix1"),
                          output = c("top1", "topk", "conformal"),
                          k = 3L) {
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

  n    <- nrow(S)
  dims <- colnames(S)
  Y    <- to_logit(S)
  idx  <- if (is.null(flagged)) seq_len(n) else which(flagged)
  n_top <- min(5L, ncol(S))
  top <- t(vapply(seq_len(n), function(i)
    dims[order(S[i, ], decreasing = TRUE)[seq_len(n_top)]],
    character(n_top)))

  pm_names <- names(pair_models)
  LL <- matrix(-Inf, n, length(pm_names), dimnames = list(NULL, pm_names))
  for (pl in pm_names) LL[, pl] <- mvn_ll_rows(Y, pair_models[[pl]])

  rank_cell <- function(i) {
    cand <- make_candidates(top[i, ])[[candidates]]
    pls  <- pair_label(cand$c1, cand$c2)
    lls  <- rep(-Inf, nrow(cand))
    in_m <- pls %in% pm_names
    if (any(in_m)) lls[in_m] <- LL[i, pls[in_m]]
    ord <- order(-lls, seq_along(lls))
    list(cand = cand, ord = ord[in_m[ord]], lls = lls, has_model = any(in_m))
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
  }
  out <- Filter(NROW, out)
  if (!length(out)) empty else do.call(rbind, out)
}

#' Identify doublet composition with an optional detector
#'
#' Scores a count matrix against one or more labelled references, learns
#' logit-normal models for synthetic heterotypic doublets, and assigns candidate
#' cell-type pairs. The built-in detector combines PhiSpace score features with
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
#' @param output One top pair, the top `k` pairs, or a conformal set.
#' @param k Number of pairs returned when `output = "topk"`.
#' @param alpha Miscoverage level for conformal sets.
#' @param seed Non-negative integer controlling simulation, scoring, and training.
#' @param score_fn Advanced scoring function with the same interface as
#'   `score_phispace()`.
#' @return An `obligato_composition` object containing composition calls,
#'   detection scores and flags, fitted composition models, parameters, and provenance.
#' @export
compose <- function(query, reference, phenotypes, hc_singlets,
                    hc_labels = NULL, detector = "builtin",
                    detector_threshold = NULL, n_per_pair = 200L,
                    candidates = c("all5", "flex2", "fix1"),
                    output = c("top1", "topk", "conformal"),
                    k = 3L, alpha = 0.10, seed = 1L,
                    score_fn = score_phispace) {
  .validate_counts(query)
  candidates <- match.arg(candidates)
  output <- match.arg(output)
  n_per_pair <- .validate_scalar_integer(n_per_pair, "n_per_pair", 30L)
  seed <- .validate_scalar_integer(seed, "seed", 0L)
  k <- .validate_scalar_integer(k, "k", 1L)
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
  if (!length(mods$pair_models)) {
    stop("No doublet-pair models could be fitted; increase `n_per_pair`.", call. = FALSE)
  }

  det_score <- threshold <- flag <- detection_model <- NULL
  if (identical(detector, "builtin")) {
    if (!length(mods$sing_models)) {
      stop("No singlet models could be fitted; provide at least 15 cells for one type.",
           call. = FALSE)
    }
    n_levels <- length(phenotypes)
    suffix <- function(i) if (n_levels > 1L) paste0("_l", i) else ""
    level_features <- function(i) {
      if (i == 1L) {
        Sa <- S_all; Ss <- S_sing; Sd <- S_dbl; m <- mods
      } else {
        Sa <- .with_seed(seed + 100L + i,
          score_fn(cnt, refs[[i]], phenotypes[[i]], sprintf(" all-cells L%d", i)))
        Sa <- .validate_scores(Sa, ncol(cnt), sprintf("Query scores for level %d", i))
        Sj <- .with_seed(seed + 200L + i,
          score_fn(joint_counts, refs[[i]], phenotypes[[i]],
                   sprintf(" train-joint L%d", i)))
        Sj <- .validate_scores(Sj, ncol(joint_counts),
                               sprintf("Training scores for level %d", i),
                               expected = colnames(Sa))
        Ss <- Sj[seq_len(n_sing), , drop = FALSE]
        Sd <- Sj[n_sing + seq_len(ncol(sim$counts)), , drop = FALSE]
        m <- fit_doublet_models(Ss, hc_labels, Sd, sim$pair_label, types)
        if (!length(m$sing_models) || !length(m$pair_models)) {
          stop(sprintf("No usable models could be fitted for annotation level %d.", i),
               call. = FALSE)
        }
      }
      list(
        sing = compute_features(Ss, m$sing_models, m$pair_models, suffix(i)),
        dbl = compute_features(Sd, m$sing_models, m$pair_models, suffix(i)),
        all = compute_features(Sa, m$sing_models, m$pair_models, suffix(i)),
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
                          output = output, k = k)
  }

  detector_name <- if (is.character(detector)) detector else "byo"
  info <- list(
    n_cells = ncol(cnt), n_hc_singlets = n_sing, types = types,
    phenotypes = phenotypes, detector = detector_name,
    detection_features = if (identical(detector, "builtin"))
      "PhiSpace scores plus library size" else NULL,
    candidates = candidates, output = output, n_per_pair = n_per_pair,
    seed = seed, threshold = threshold,
    alpha = if (output == "conformal") alpha else NULL)

  structure(
    list(composition = comp, detection_score = det_score, flag = flag,
         conformal = conformal, models = mods,
         detection_model = detection_model, info = info,
         provenance = .package_provenance()),
    class = "obligato_composition")
}
