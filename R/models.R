##############################################################################
# Logit-normal pair models.
# Extracted VERBATIM from run_methodC_compose.R (M1). Do not modify behaviour here —
# see package_plan.md §10: refactor and behaviour change must never share a commit.
##############################################################################

to_logit <- function(x, eps = 1e-3) {
  u <- (pmin(pmax(x, -1 + eps), 1 - eps) + 1) / 2
  log(u / (1 - u))
}

pair_label <- function(a, b) ifelse(a <= b, paste(a, "+", b), paste(b, "+", a))

fit_ln_shrink <- function(Y) {
  mu <- colMeans(Y)
  Sigma <- as.matrix(suppressMessages(cov.shrink(Y, verbose = FALSE)))
  ch <- chol(Sigma)
  list(mu = mu, Sigma_inv = chol2inv(ch), log_det = 2 * sum(log(diag(ch))))
}

quad_form_rows <- function(Y, mu, Sigma_inv) {
  D <- sweep(Y, 2, mu, "-")
  rowSums((D %*% Sigma_inv) * D)
}

mvn_ll_rows <- function(Y, model) {
  p <- ncol(Y)
  -0.5 * (p * log(2 * pi) + model$log_det + quad_form_rows(Y, model$mu, model$Sigma_inv))
}

marginal_ll <- function(y, mu, sigma) {
  -0.5 * log(2 * pi) - log(pmax(sigma, 1e-6)) -
    0.5 * ((y - mu) / pmax(sigma, 1e-6))^2
}

##############################################################################
# Fit the per-type singlet models and per-pair doublet models (M1b).
#
# Extracted from run_compose_kingetal.R:227-247 (Section 3), a block that is
# byte-identical to run_methodC_compose.R's process_batch() fitting loop, so one
# function serves both the King case study and the benchmark (needed for the M2b
# replication gate).
#
# Fits, in logit space:
#   * one shrunk logit-normal per singlet type with >= `min_singlet` cells;
#   * per-type marginal SDs (consumed by marginal_ll);
#   * one shrunk logit-normal per simulated pair with >= `min_pair` doublets.
#
# Behaviour-preserving deviations from the inline block, all documented:
#   * the hardcoded 15 / 30 minimum-cell cutoffs become `min_singlet` / `min_pair`
#     arguments defaulting to the same values;
#   * the script's ambient `types` vector becomes the `types` argument, defaulting
#     to the sorted unique singlet types. Iteration order affects only the order of
#     the returned named lists, which is numerically irrelevant downstream
#     (compute_features and mvn_ll_rows reduce over models with max);
#   * the script's cat() progress logging is dropped.
# For identical inputs (and matching `types`) every fitted model is identical to
# the script's.
fit_doublet_models <- function(S_train_sing, hc_types, S_train_dbl, dbl_pl,
                               types = sort(unique(hc_types)),
                               min_singlet = 15, min_pair = 30) {
  Y_ts <- to_logit(S_train_sing)
  sing_models <- list()
  for (tt in types) {
    idx <- which(hc_types == tt)
    if (length(idx) >= min_singlet) {
      sing_models[[tt]] <- fit_ln_shrink(Y_ts[idx, , drop = FALSE])
    }
  }

  sing_marginal_sd <- list()
  for (tt in names(sing_models)) {
    idx <- which(hc_types == tt)
    sing_marginal_sd[[tt]] <- apply(Y_ts[idx, , drop = FALSE], 2, sd)
  }

  Y_td <- to_logit(S_train_dbl)
  pair_models <- list()
  for (pl in unique(dbl_pl)) {
    idx <- which(dbl_pl == pl)
    if (length(idx) >= min_pair) {
      pair_models[[pl]] <- fit_ln_shrink(Y_td[idx, , drop = FALSE])
    }
  }

  list(sing_models = sing_models,
       sing_marginal_sd = sing_marginal_sd,
       pair_models = pair_models)
}
