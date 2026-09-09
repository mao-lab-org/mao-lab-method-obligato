# Internal split-conformal prediction sets for composition.
#
# Negative pair log-likelihood is the nonconformity score. Under exchangeability,
# the finite-sample order-statistic threshold covers the true pair with probability
# at least 1 - alpha. compose() uses an independent synthetic calibration sample;
# empirical coverage should be checked whenever labelled doublets are available.--- .gitignore

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

  sorted <- sort(s_cal)
  q <- vapply(alpha, function(a) {
    rank <- min(ceiling((1 - a) * (n_cal + 1)), n_cal)
    sorted[[rank]]
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
