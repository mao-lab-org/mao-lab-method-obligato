# High-confidence singlet selection with a numerically guarded logit transform.
.EPS <- 1e-7

.clip01 <- function(s) pmin(pmax(s, .EPS), 1 - .EPS)

.logit  <- function(s) { sc <- .clip01(s); log(sc / (1 - sc)) }

fit_logitnorm_mix <- function(s, init_class, n_iter = 200L, tol = 1e-7) {
  if (!is.numeric(s) || !length(s) || length(init_class) != length(s) ||
      anyNA(s) || any(!is.finite(s)) || any(s < 0 | s > 1) || anyNA(init_class) ||
      any(!init_class %in% 1:2) || any(tabulate(init_class, 2L) < 2L)) {
    stop("`s` must be finite; `init_class` must assign at least two values to each class.",
         call. = FALSE)
  }
  n_iter <- .validate_scalar_integer(n_iter, "n_iter", 1L)
  if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive number.", call. = FALSE)
  }

  y <- .logit(s); N <- length(y); K <- 2L
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
      sigma[k] <- max(sqrt(sum(w * (y - mu[k])^2) / W), 1e-6)
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

#' Estimate a high-confidence singlet cutoff
#'
#' Fits a two-component logit-normal mixture to a score and returns the score at
#' which the posterior probability of the higher-score component reaches the
#' requested value.
#' @param s Numeric scores in the unit interval.
#' @param init_class Initial component assignments coded as 1 and 2, with at least
#'   two observations in each component.
#' @param target_posterior Posterior probability used to define the cutoff.
#' @return A list containing `cutoff` and the fitted mixture `fit`.
#' @export
hc_singlet_cutoff <- function(s, init_class, target_posterior = 0.01) {
  if (length(target_posterior) != 1L || !is.numeric(target_posterior) ||
      is.na(target_posterior) || target_posterior <= 0 || target_posterior >= 1) {
    stop("`target_posterior` must be strictly between zero and one.", call. = FALSE)
  }
  fit <- fit_logitnorm_mix(s, init_class)
  f <- function(x) posterior_doublet_score(x, fit) - target_posterior
  if (f(.EPS) > 0) return(list(cutoff = .EPS, fit = fit))
  if (f(1 - .EPS) < 0) return(list(cutoff = 1 - .EPS, fit = fit))
  s_star <- uniroot(f, interval = c(.EPS, 1 - .EPS), tol = 1e-8)$root
  list(cutoff = s_star, fit = fit)
}
