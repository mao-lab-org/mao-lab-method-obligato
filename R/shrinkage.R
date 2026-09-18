# Default empirical-Bayes weight for the detection features. Fixed rather than
# cross-validated: the detection response is flat over roughly 0.7-0.9, and 0.75
# keeps a margin from lambda = 1, where the pair's own simulated doublets stop
# contributing and within-lineage detection degrades sharply.
.EB_LAMBDA_DEFAULT <- 0.75

# Empirical-Bayes shrinkage of the doublet-pair covariances (M2f).
#
# Each pair covariance is blended toward a structured target built from the two
# constituents' SINGLET models:
#
#   Sigma_eb = (1 - lambda) * Sigma_lw + lambda * Sigma_0.
#
# The target follows from how a doublet is formed. Its counts are c_a + c_b, so
# after depth normalisation its score vector is approximately a random convex
# combination pi * x_a + (1 - pi) * x_b, with pi the depth share of cell a. That
# mixture is linear in RAW score space, whereas the models live in logit space,
# so the target is obtained by simulating the mixture and transforming, which
# carries the nonlinearity exactly. In raw space the mixture has closed-form
# moments (see .mixture_moments_raw), which the simulation is tested against.
#
# lambda is chosen by K-fold exact-pair accuracy on the SIMULATED training
# doublets, whose pairs are known by construction, so no ground truth is needed
# at run time. Selecting it by held-out likelihood instead optimises the wrong
# objective: it drives lambda toward 1, where composition accuracy degrades.
#
# Scope, deliberately fixed: singlet models are NOT shrunk (they carry the
# reference-only log_det diagnostic), and lambda = 0 is always in the grid, so
# the procedure can decline to shrink on data where shrinkage would hurt.

# Mean and Ledoit-Wolf covariance of a set of cells in RAW score space.
.raw_moments <- function(S) {
  list(mu = colMeans(S),
       Sigma = as.matrix(suppressMessages(corpcor::cov.shrink(S, verbose = FALSE))))
}

# Closed-form moments of pi * x_a + (1 - pi) * x_b in raw space, with pi
# independent of x. Used to validate the Monte-Carlo target.
.mixture_moments_raw <- function(ma, Sa, mb, Sb, mu_pi, var_pi) {
  d <- ma - mb
  list(mu = mu_pi * ma + (1 - mu_pi) * mb,
       Sigma = (mu_pi^2 + var_pi) * Sa + ((1 - mu_pi)^2 + var_pi) * Sb +
         var_pi * tcrossprod(d))
}

# Empirical distribution of the depth share pi = L_a / (L_a + L_b).
.pi_pool <- function(lib_a, lib_b, n = 5000L) {
  la <- sample(lib_a, n, replace = TRUE)
  lb <- sample(lib_b, n, replace = TRUE)
  la / (la + lb)
}

.mvn_draw <- function(n, mu, Sigma) {
  d <- length(mu)
  e <- eigen(Sigma, symmetric = TRUE)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 1e-10)), d, d)
  sweep(matrix(stats::rnorm(n * d), n, d) %*% t(A), 2, mu, "+")
}

# One Monte-Carlo sample from the mixture, in RAW score space.
.mc_pair_sample <- function(ma, Sa, mb, Sb, pi_pool, R = 2000L) {
  p <- sample(pi_pool, R, replace = TRUE)
  p * .mvn_draw(R, ma, Sa) + (1 - p) * .mvn_draw(R, mb, Sb)
}

# The structured target: the mixture pushed through to_logit.
.mc_pair_target <- function(ma, Sa, mb, Sb, pi_pool, R = 2000L) {
  Z <- to_logit(.mc_pair_sample(ma, Sa, mb, Sb, pi_pool, R))
  list(mu = colMeans(Z), Sigma = stats::cov(Z))
}

# Convex combination, returned in the same shape fit_ln_shrink produces.
.blend_pair_model <- function(Y, target, lambda) {
  if (lambda <= 0) return(NULL)
  S_lw <- as.matrix(suppressMessages(corpcor::cov.shrink(Y, verbose = FALSE)))
  Sigma <- (1 - lambda) * S_lw + lambda * target$Sigma
  mu <- (1 - lambda) * colMeans(Y) + lambda * target$mu
  ch <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  list(mu = mu, Sigma_inv = chol2inv(ch), log_det = 2 * sum(log(diag(ch))))
}

# Targets for every pair that has both constituents' singlet models.
.pair_targets <- function(S_sing, sing_labels, lib_sing, pair_names,
                          R = 2000L, min_singlet = 15L) {
  types <- sort(unique(sing_labels))
  raw <- list(); pool <- list()
  for (tt in types) {
    idx <- which(sing_labels == tt)
    if (length(idx) < min_singlet) next
    raw[[tt]] <- .raw_moments(S_sing[idx, , drop = FALSE])
    pool[[tt]] <- lib_sing[idx]
  }
  out <- list()
  for (pl in pair_names) {
    ab <- strsplit(pl, " + ", fixed = TRUE)[[1L]]
    if (length(ab) != 2L || is.null(raw[[ab[1L]]]) || is.null(raw[[ab[2L]]])) next
    out[[pl]] <- .mc_pair_target(raw[[ab[1L]]]$mu, raw[[ab[1L]]]$Sigma,
                                 raw[[ab[2L]]]$mu, raw[[ab[2L]]]$Sigma,
                                 .pi_pool(pool[[ab[1L]]], pool[[ab[2L]]]), R)
  }
  out
}

# K-fold exact-pair accuracy on the simulated doublets, for one lambda. Folds
# are assigned WITHIN each pair so every pair model is still fitted, and the
# held-out doublets of all pairs are pooled for scoring.
.cv_lambda_accuracy <- function(S_dbl, dbl_pair, targets, lambda, fold,
                                min_pair = 30L) {
  pls <- intersect(unique(dbl_pair), names(targets))
  accs <- vapply(sort(unique(fold)), function(f) {
    pm <- list()
    for (pl in pls) {
      idx <- which(dbl_pair == pl & fold != f)
      if (length(idx) < min_pair) next
      Y <- to_logit(S_dbl[idx, , drop = FALSE])
      m <- if (lambda <= 0) fit_ln_shrink(Y) else
        .blend_pair_model(Y, targets[[pl]], lambda)
      if (!is.null(m)) pm[[pl]] <- m
    }
    te <- which(fold == f & dbl_pair %in% names(pm))
    if (!length(te) || !length(pm)) return(NA_real_)
    # Scored with the same candidate rule compose() uses, so the selected weight
    # is chosen under the criterion that will actually be applied. compose_pairs()
    # errors if the model names and score columns disagree; it is not this
    # function's job to paper over that.
    cp <- compose_pairs(S_dbl[te, , drop = FALSE], pm,
                        candidates = "all5", output = "top1")
    pred <- pair_label(cp$c1, cp$c2)
    mean(pred == dbl_pair[te], na.rm = TRUE)
  }, numeric(1))
  mean(accs, na.rm = TRUE)
}

#' Empirical-Bayes shrinkage of fitted pair models
#'
#' Blends each pair covariance toward a target predicted from the two
#' constituents' singlet models. Singlet models are returned unchanged.
#'
#' @param models Output of `fit_doublet_models()`.
#' @param S_sing,sing_labels Training singlet scores and their type labels.
#' @param S_dbl,dbl_pair Simulated training doublet scores and their pair labels.
#' @param lib_sing Library sizes of the training singlets, for the depth-share
#'   distribution.
#' @param lambda `"cv"` to select by cross-validation, or a number in `[0, 1]`.
#' @param grid Candidate values when `lambda = "cv"`. Must contain 0.
#' @param R Monte-Carlo draws per target.
#' @param nfold Cross-validation folds.
#' @param seed Seed for folds and Monte-Carlo draws.
#' @return A list with the updated `models`, the chosen `lambda`, and the `cv`
#'   curve (`NULL` when lambda was supplied).
#' @keywords internal
apply_eb_pairs <- function(models, S_sing, sing_labels, S_dbl, dbl_pair,
                           lib_sing, lambda = "cv",
                           grid = c(0, 0.25, 0.5, 0.75), R = 2000L,
                           nfold = 5L, seed = 1L, min_pair = 30L) {
  stopifnot(is.list(models), !is.null(models$pair_models))
  if (!identical(lambda, "cv")) {
    if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
        lambda < 0 || lambda > 1) {
      stop("`lambda` must be \"cv\" or a single number in [0, 1].", call. = FALSE)
    }
  }
  if (!0 %in% grid) {
    stop("The lambda grid must contain 0, so that shrinkage can be declined.",
         call. = FALSE)
  }
  pair_names <- names(models$pair_models)
  if (!length(pair_names)) return(list(models = models, lambda = 0, cv = NULL))

  targets <- .with_seed(seed, .pair_targets(S_sing, sing_labels, lib_sing,
                                            pair_names, R = R))
  if (!length(targets)) return(list(models = models, lambda = 0, cv = NULL))

  cv <- NULL
  if (identical(lambda, "cv")) {
    fold <- .with_seed(seed + 1L, {
      f <- integer(length(dbl_pair))
      for (pl in unique(dbl_pair)) {
        i <- which(dbl_pair == pl)
        f[i] <- sample(rep_len(seq_len(nfold), length(i)))
      }
      f
    })
    cv <- vapply(grid, function(l)
      .cv_lambda_accuracy(S_dbl, dbl_pair, targets, l, fold, min_pair),
      numeric(1))
    names(cv) <- as.character(grid)
    # Too little data for any fold to be scored: decline to shrink rather than
    # picking arbitrarily. This is the same conservative default as a CV that
    # genuinely prefers lambda = 0.
    lambda <- if (all(!is.finite(cv))) 0 else
      grid[which.max(replace(cv, !is.finite(cv), -Inf))]
  }
  if (lambda <= 0) return(list(models = models, lambda = 0, cv = cv))

  n_done <- 0L
  for (pl in pair_names) {
    if (is.null(targets[[pl]])) next
    idx <- which(dbl_pair == pl)
    if (length(idx) < min_pair) next
    m <- .blend_pair_model(to_logit(S_dbl[idx, , drop = FALSE]),
                           targets[[pl]], lambda)
    if (!is.null(m)) { models$pair_models[[pl]] <- m; n_done <- n_done + 1L }
  }
  list(models = models, lambda = lambda, cv = cv, n_shrunk = n_done)
}

#' Cross-validate the EB weight on a detection objective
#'
#' The composition criterion in [.cv_lambda_accuracy()] selects `lambda` on
#' exact-pair accuracy. When the shrunk models are used only to build detection
#' features -- which is what [compose()] does -- that optimises a task the weight
#' is not applied to. This is the matching criterion: each held-out fold is scored
#' with detection AUPRC, using the same feature block and classifier the shipped
#' pipeline uses.
#'
#' Pair models are refitted inside every fold, on the training folds only.
#' Selecting `lambda` against models that had already seen the held-out simulated
#' doublets would leak, and would systematically flatter larger `lambda`.
#'
#' Be aware of what this criterion can and cannot resolve. Measured on the three
#' benchmark datasets, it detects reliably that shrinkage helps -- the deficit at
#' `lambda = 0` is around four times the fold-to-fold standard deviation -- but
#' across `lambda >= 0.5` the spread is roughly the size of that standard
#' deviation, so the arg-max in that region is close to arbitrary. It is offered
#' as an option, not as the default.
#'
#' @param levels List with one entry per annotation level, each a list with
#'   `Ss` (singlet scores) and `Sd` (simulated-doublet scores).
#' @param sing_labels,dbl_pair Labels for the singlets and the simulated pairs.
#' @param types Reference type vocabulary.
#' @param cnt_sing,cnt_dbl Counts for the singlets and simulated doublets, for
#'   the library-size feature block.
#' @param lib_sing Singlet library sizes, for the EB target.
#' @param grid Candidate weights; must contain 0.
#' @param R,nfold,seed,min_pair As in [apply_eb_pairs()].
#' @param suffix Function mapping a level index to its feature-name suffix.
#' @return A list with the selected `lambda` and the `cv` curve over `grid`.
#' @keywords internal
#' @noRd
.cv_lambda_detection <- function(levels, sing_labels, dbl_pair, types,
                                 cnt_sing, cnt_dbl, lib_sing, grid,
                                 R = 2000L, nfold = 5L, seed = 1L,
                                 suffix = function(i) "", min_pair = 30L) {
  n_s <- nrow(levels[[1L]]$Ss); n_d <- nrow(levels[[1L]]$Sd)
  fold_s <- .with_seed(seed + 2L, sample(rep_len(seq_len(nfold), n_s)))
  fold_d <- .with_seed(seed + 3L, sample(rep_len(seq_len(nfold), n_d)))

  one_fold <- function(lambda, f) {
    tr_s <- which(fold_s != f); te_s <- which(fold_s == f)
    tr_d <- which(fold_d != f); te_d <- which(fold_d == f)
    if (!length(te_s) || !length(te_d)) return(NA_real_)
    F <- list(str = list(), dtr = list(), ste = list(), dte = list())
    for (i in seq_along(levels)) {
      L <- levels[[i]]
      m <- fit_doublet_models(L$Ss[tr_s, , drop = FALSE], sing_labels[tr_s],
                              L$Sd[tr_d, , drop = FALSE], dbl_pair[tr_d], types)
      if (!length(m$sing_models) || !length(m$pair_models)) return(NA_real_)
      if (lambda > 0) {
        m <- apply_eb_pairs(m, L$Ss[tr_s, , drop = FALSE], sing_labels[tr_s],
                            L$Sd[tr_d, , drop = FALSE], dbl_pair[tr_d],
                            lib_sing[tr_s], lambda = lambda, grid = grid,
                            R = R, seed = seed, min_pair = min_pair)$models
      }
      sfx <- suffix(i)
      F$str[[i]] <- compute_features(L$Ss[tr_s, , drop = FALSE], m$sing_models, m$pair_models, sfx)
      F$dtr[[i]] <- compute_features(L$Sd[tr_d, , drop = FALSE], m$sing_models, m$pair_models, sfx)
      F$ste[[i]] <- compute_features(L$Ss[te_s, , drop = FALSE], m$sing_models, m$pair_models, sfx)
      F$dte[[i]] <- compute_features(L$Sd[te_d, , drop = FALSE], m$sing_models, m$pair_models, sfx)
    }
    S1 <- levels[[1L]]$Ss; D1 <- levels[[1L]]$Sd
    # the type-median baseline is a training-fold quantity too
    tm <- libsize_type_medians(cnt_sing[, tr_s, drop = FALSE], S1[tr_s, , drop = FALSE])
    lb <- function(x, sc) libsize_features(x, top1_type(sc), tm)
    f_str <- cbind(do.call(cbind, F$str), lb(cnt_sing[, tr_s, drop = FALSE], S1[tr_s, , drop = FALSE]))
    f_dtr <- cbind(do.call(cbind, F$dtr), lb(cnt_dbl[, tr_d, drop = FALSE], D1[tr_d, , drop = FALSE]))
    f_ste <- cbind(do.call(cbind, F$ste), lb(cnt_sing[, te_s, drop = FALSE], S1[te_s, , drop = FALSE]))
    f_dte <- cbind(do.call(cbind, F$dte), lb(cnt_dbl[, te_d, drop = FALSE], D1[te_d, , drop = FALSE]))
    bst <- train_detector(f_str, f_dtr, seed = seed)
    p <- as.numeric(stats::predict(bst, as.matrix(rbind(f_ste, f_dte))))
    y <- c(rep(0L, nrow(f_ste)), rep(1L, nrow(f_dte)))
    compute_auprc(y, p)$auprc
  }

  cv <- vapply(grid, function(l)
    mean(vapply(seq_len(nfold), function(f) one_fold(l, f), numeric(1)), na.rm = TRUE),
    numeric(1))
  names(cv) <- as.character(grid)
  lambda <- if (all(!is.finite(cv))) 0 else grid[which.max(replace(cv, !is.finite(cv), -Inf))]
  list(lambda = lambda, cv = cv)
}
