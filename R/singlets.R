##############################################################################
# High-confidence singlet selection.
# Extracted VERBATIM from eval_harness.R (M1).
##############################################################################

# Numerical guard for the logit transform. Extracted verbatim from
# eval_harness.R:56 — omitted in the first M1 pass, which the parity gate caught
# (`object '.EPS' not found`). Exactly the class of hidden dependency the gate exists for.
.EPS <- 1e-7

.clip01 <- function(s) pmin(pmax(s, .EPS), 1 - .EPS)

.logit  <- function(s) { sc <- .clip01(s); log(sc / (1 - sc)) }

fit_logitnorm_mix <- function(s, init_class, n_iter = 200, tol = 1e-7) {

  y <- .logit(s); N <- length(y); K <- 2
  mu <- sigma <- pi_k <- numeric(K)
  for (k in 1:K) {
    yi <- y[init_class == k]
    mu[k]    <- mean(yi)
    sigma[k] <- max(sd(yi), 1e-3)
    pi_k[k]  <- length(yi) / N
  }

  ll_prev <- -Inf
  for (iter in seq_len(n_iter)) {
    log_lik_k <- cbind(
      log(pi_k[1]) + dnorm(y, mu[1], sigma[1], log = TRUE),
      log(pi_k[2]) + dnorm(y, mu[2], sigma[2], log = TRUE)
    )
    log_max  <- pmax(log_lik_k[, 1], log_lik_k[, 2])
    log_norm <- log_max + log(rowSums(exp(log_lik_k - log_max)))
    gamma    <- exp(log_lik_k - log_norm)

    pi_k <- colMeans(gamma)
    for (k in 1:K) {
      w  <- gamma[, k]
      W  <- sum(w)
      mu[k]    <- sum(w * y) / W
      sigma[k] <- sqrt(sum(w * (y - mu[k])^2) / W)
    }

    ll <- sum(log_norm)
    if (abs(ll - ll_prev) < tol) break
    ll_prev <- ll
  }

  if (mu[1] > mu[2]) {
    mu <- mu[2:1]; sigma <- sigma[2:1]; pi_k <- pi_k[2:1]
  }
  list(mu = mu, sigma = sigma, pi = pi_k, logLik = ll)
}

posterior_doublet_score <- function(s, fit) {
  y  <- .logit(s)
  l1 <- log(fit$pi[1]) + dnorm(y, fit$mu[1], fit$sigma[1], log = TRUE)
  l2 <- log(fit$pi[2]) + dnorm(y, fit$mu[2], fit$sigma[2], log = TRUE)
  m  <- pmax(l1, l2)
  exp(l2 - m) / (exp(l1 - m) + exp(l2 - m))
}

# Cutoff on raw score s such that P(doublet | s*) = target.
hc_singlet_cutoff <- function(s, init_class, target_posterior = 0.01) {
  fit <- fit_logitnorm_mix(s, init_class)
  f <- function(x) posterior_doublet_score(x, fit) - target_posterior
  if (f(.EPS) > 0) return(list(cutoff = .EPS, fit = fit))
  if (f(1 - .EPS) < 0) return(list(cutoff = 1 - .EPS, fit = fit))
  s_star <- uniroot(f, interval = c(.EPS, 1 - .EPS), tol = 1e-8)$root
  list(cutoff = s_star, fit = fit)
}
