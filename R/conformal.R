##############################################################################
# Split-conformal prediction sets for doublet composition (S4 / Item 3).
#
# Extracted from Phase0_Detection/conformal_composition.R (Sections 3-5), the
# base-R split-conformal procedure with negative-log-likelihood nonconformity.
# Factored into two pure functions so the algorithm is golden-testable in
# isolation and compose() can wire it to whatever calibration source is available.
#
# Nonconformity of a candidate pair p for a cell is s = -loglik_p(cell): the true
# pair should score high (small s). Given a labelled calibration set we take the
# (1-alpha)-quantile of the true-pair nonconformity as q_hat, and the prediction
# set for a query cell is every pair whose nonconformity <= q_hat, i.e. every pair
# with loglik >= -q_hat. With exchangeable calibration/query draws this set covers
# the true pair with probability >= 1 - alpha.
#
# NOTE ON CALIBRATION SOURCE. The eval script calibrates on ground-truth doublets
# (known pairs) held out from the test set. A production query has no ground-truth
# pairs, so compose() instead calibrates on freshly simulated doublets, whose pairs
# are known by construction and are drawn from the same per-pair models. The pure
# functions below are agnostic to which source supplies (log_h, true_pairs).
##############################################################################

# All-pairs log-likelihood matrix: cell x pair-model, colnames = pair labels.
# Shared by compose_pairs' ranking and the conformal path (mvn_ll_rows batched).
all_pair_ll <- function(S, pair_models) {
  Y  <- to_logit(S)
  pm <- names(pair_models)
  LL <- matrix(-Inf, nrow(S), length(pm), dimnames = list(NULL, pm))
  for (pl in pm) LL[, pl] <- mvn_ll_rows(Y, pair_models[[pl]])
  LL
}

# Calibrate the conformal threshold(s) q_hat from a labelled calibration set.
#   log_h      : n_cal x P matrix of pair log-likelihoods (colnames = pair labels)
#   true_pairs : length-n_cal character vector of the true pair per calibration cell
#   alpha      : miscoverage level(s); 1 - alpha is the target coverage
# Returns a named numeric vector of q_hat, one per alpha. A true pair absent from
# the fitted models contributes nonconformity Inf (it can never be covered), which
# is the honest, conservative behaviour.
conformal_calibrate <- function(log_h, true_pairs, alpha = 0.10) {
  stopifnot(is.matrix(log_h), nrow(log_h) == length(true_pairs),
            all(alpha > 0), all(alpha < 1))
  pm    <- colnames(log_h)
  n_cal <- nrow(log_h)
  if (n_cal == 0L) stop("empty calibration set", call. = FALSE)

  s_cal <- vapply(seq_len(n_cal), function(i) {
    tp <- true_pairs[i]
    if (!is.na(tp) && tp %in% pm) -log_h[i, tp] else Inf
  }, numeric(1))

  q <- vapply(alpha, function(a) {
    level <- min(ceiling((1 - a) * (n_cal + 1)) / n_cal, 1)  # split-conformal level
    unname(quantile(s_cal, probs = level, type = 7))
  }, numeric(1))
  names(q) <- paste0("alpha_", alpha)
  q
}

# Build prediction sets for query cells at a single threshold q_hat.
#   log_h : n_query x P matrix of pair log-likelihoods (colnames = pair labels)
#   q_hat : scalar nonconformity threshold from conformal_calibrate
# Returns a list with the per-cell pair vector (`sets`) and integer `sizes`.
conformal_sets <- function(log_h, q_hat) {
  stopifnot(is.matrix(log_h), length(q_hat) == 1L)
  pm     <- colnames(log_h)
  in_set <- log_h >= -q_hat                       # membership: nonconformity <= q_hat
  sets   <- lapply(seq_len(nrow(log_h)), function(i) pm[which(in_set[i, ])])
  list(sets = sets, sizes = lengths(sets))
}

# Orchestrator: calibrate on (log_h_cal, true_pairs_cal), then set-build for the
# query rows of log_h_query, at one alpha. Returns sets, sizes, and q_hat.
conformal_pairs <- function(log_h_cal, true_pairs_cal, log_h_query, alpha = 0.10) {
  stopifnot(length(alpha) == 1L)
  q_hat <- conformal_calibrate(log_h_cal, true_pairs_cal, alpha = alpha)[[1]]
  cs    <- conformal_sets(log_h_query, q_hat)
  list(sets = cs$sets, sizes = cs$sizes, q_hat = q_hat, alpha = alpha)
}
