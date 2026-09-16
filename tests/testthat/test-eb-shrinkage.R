## M2f — empirical-Bayes pair covariances. The tests are organised around the
## claims the implementation has to honour, not around the functions.

set.seed(11)

# Small synthetic world. The score columns ARE the reference types, as they are
# in real data -- `compose_pairs(candidates = "all5")` builds its candidate pairs
# from the score column names, so a fixture that names them otherwise silently
# matches nothing.
types <- c("A", "B", "C", "D", "E", "F")
d <- length(types)
mk_sing <- function(n, centre) {
  m <- matrix(stats::rnorm(n * d, 0, 0.25), n, d)
  m <- sweep(m, 2, centre, "+")
  matrix(pmin(pmax(m, -0.98), 0.98), n, d, dimnames = list(NULL, types))
}
# Each type scores high on its own column.
centres <- lapply(seq_along(types), function(i) {
  v <- rep(-0.1, d); v[i] <- 0.75; v
})
names(centres) <- types
n_per <- 120L
S_sing <- do.call(rbind, lapply(types, function(tt) mk_sing(n_per, centres[[tt]])))
sing_labels <- rep(types, each = n_per)
lib_sing <- stats::rlnorm(nrow(S_sing), log(3000), 0.3)

# Simulated doublets: mix two singlets in raw space, as the real generator does
# after normalisation, then clip back into range.
mk_dbl <- function(a, b, n) {
  ia <- sample(which(sing_labels == a), n, TRUE)
  ib <- sample(which(sing_labels == b), n, TRUE)
  p <- stats::runif(n, 0.35, 0.65)
  m <- p * S_sing[ia, , drop = FALSE] + (1 - p) * S_sing[ib, , drop = FALSE]
  matrix(pmin(pmax(m, -0.98), 0.98), n, d, dimnames = list(NULL, colnames(S_sing)))
}
pairs_all <- utils::combn(types, 2L, simplify = FALSE)
n_dbl <- 90L
S_dbl <- do.call(rbind, lapply(pairs_all, function(p) mk_dbl(p[1], p[2], n_dbl)))
dbl_pair <- unlist(lapply(pairs_all, function(p)
  rep(Obligato:::pair_label(p[1], p[2]), n_dbl)))

base <- Obligato:::fit_doublet_models(S_sing, sing_labels, S_dbl, dbl_pair, types)

recover_Sigma <- function(m) solve(m$Sigma_inv)

test_that("the Monte-Carlo target matches the closed-form mixture moments", {
  # The target is built by simulation because the logit transform is nonlinear,
  # but in RAW space the mixture has closed-form moments. This checks the
  # simulation against that algebra, which is the only independent check of the
  # target's correctness.
  ma <- Obligato:::.raw_moments(S_sing[sing_labels == "A", , drop = FALSE])
  mb <- Obligato:::.raw_moments(S_sing[sing_labels == "B", , drop = FALSE])
  pi_pool <- stats::runif(20000, 0.2, 0.8)
  set.seed(4)
  X <- Obligato:::.mc_pair_sample(ma$mu, ma$Sigma, mb$mu, mb$Sigma,
                                  pi_pool, R = 200000L)
  cf <- Obligato:::.mixture_moments_raw(ma$mu, ma$Sigma, mb$mu, mb$Sigma,
                                        mean(pi_pool), stats::var(pi_pool))
  expect_equal(colMeans(X), cf$mu, tolerance = 0.02, ignore_attr = TRUE)
  expect_equal(stats::cov(X), cf$Sigma, tolerance = 0.05, ignore_attr = TRUE)
})

test_that("a degenerate depth share reduces the mixture to one constituent", {
  ma <- Obligato:::.raw_moments(S_sing[sing_labels == "A", , drop = FALSE])
  mb <- Obligato:::.raw_moments(S_sing[sing_labels == "B", , drop = FALSE])
  set.seed(5)
  X <- Obligato:::.mc_pair_sample(ma$mu, ma$Sigma, mb$mu, mb$Sigma,
                                  rep(1, 100), R = 50000L)
  expect_equal(colMeans(X), ma$mu, tolerance = 0.02, ignore_attr = TRUE)
  expect_equal(stats::cov(X), ma$Sigma, tolerance = 0.05, ignore_attr = TRUE)
})

test_that("lambda = 0 leaves the models exactly unchanged", {
  out <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                   lib_sing, lambda = 0, R = 200L, seed = 3L)
  expect_identical(out$models, base)
  expect_identical(out$lambda, 0)
})

test_that("lambda = 1 replaces the covariance with the target exactly", {
  set.seed(7)
  tg <- Obligato:::.with_seed(3L, Obligato:::.pair_targets(
    S_sing, sing_labels, lib_sing, names(base$pair_models), R = 400L))
  out <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                   lib_sing, lambda = 1, R = 400L, seed = 3L)
  pl <- names(base$pair_models)[1]
  expect_equal(recover_Sigma(out$models$pair_models[[pl]]), unname(tg[[pl]]$Sigma),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(unname(out$models$pair_models[[pl]]$mu), unname(tg[[pl]]$mu),
               tolerance = 1e-8)
})

test_that("intermediate lambda is the stated convex combination", {
  lam <- 0.4
  tg <- Obligato:::.with_seed(3L, Obligato:::.pair_targets(
    S_sing, sing_labels, lib_sing, names(base$pair_models), R = 400L))
  out <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                   lib_sing, lambda = lam, R = 400L, seed = 3L)
  pl <- names(base$pair_models)[2]
  idx <- which(dbl_pair == pl)
  Y <- Obligato:::to_logit(S_dbl[idx, , drop = FALSE])
  S_lw <- as.matrix(suppressMessages(corpcor::cov.shrink(Y, verbose = FALSE)))
  expect_equal(recover_Sigma(out$models$pair_models[[pl]]),
               unname((1 - lam) * S_lw + lam * tg[[pl]]$Sigma),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(unname(out$models$pair_models[[pl]]$mu),
               unname((1 - lam) * colMeans(Y) + lam * tg[[pl]]$mu), tolerance = 1e-8)
})

test_that("log_det is consistent with the returned precision matrix", {
  out <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                   lib_sing, lambda = 0.5, R = 300L, seed = 3L)
  for (pl in names(out$models$pair_models)[1:3]) {
    m <- out$models$pair_models[[pl]]
    expect_equal(m$log_det, determinant(recover_Sigma(m), logarithm = TRUE)$modulus[[1]],
                 tolerance = 1e-6)
  }
})

test_that("singlet models are never touched", {
  out <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                   lib_sing, lambda = 0.75, R = 300L, seed = 3L)
  expect_identical(out$models$sing_models, base$sing_models)
  expect_identical(out$models$sing_marginal_sd, base$sing_marginal_sd)
})

test_that("selection is deterministic given a seed, and reproducible", {
  a <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                 lib_sing, lambda = "cv", grid = c(0, 0.5),
                                 R = 300L, seed = 9L)
  b <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                 lib_sing, lambda = "cv", grid = c(0, 0.5),
                                 R = 300L, seed = 9L)
  expect_identical(a$lambda, b$lambda)
  expect_equal(a$cv, b$cv)
  expect_equal(a$models$pair_models, b$models$pair_models)
})

test_that("cross-validation can and does decline to shrink when the target is wrong", {
  # Give every pair a deliberately absurd target: shrinking toward it must make
  # composition worse, so the criterion should choose lambda = 0. This is the
  # escape hatch that protected PBMC in the benchmark.
  bad <- lapply(names(base$pair_models), function(pl) {
    list(mu = rep(50, d), Sigma = diag(1e-3, d))
  })
  names(bad) <- names(base$pair_models)
  fold <- rep_len(1:3, length(dbl_pair))
  acc0 <- Obligato:::.cv_lambda_accuracy(S_dbl, dbl_pair, bad, 0, fold)
  acc_hi <- Obligato:::.cv_lambda_accuracy(S_dbl, dbl_pair, bad, 0.75, fold)
  expect_gt(acc0, acc_hi)
})

test_that("the lambda grid must offer zero", {
  expect_error(
    Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                              lib_sing, lambda = "cv", grid = c(0.25, 0.5),
                              R = 200L, seed = 3L),
    "must contain 0")
})

test_that("invalid lambda is rejected", {
  for (bad_lambda in list(-0.1, 1.5, NA_real_, "0.5", c(0.2, 0.3))) {
    expect_error(
      Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                lib_sing, lambda = bad_lambda, R = 100L),
      "lambda")
  }
})

test_that("pairs without a usable target keep their Ledoit-Wolf model", {
  # Drop one type's singlets below min_singlet: every pair containing it loses
  # its target and must fall back unchanged.
  keep <- which(sing_labels != "D" | seq_along(sing_labels) %in%
                  head(which(sing_labels == "D"), 5))
  out <- Obligato:::apply_eb_pairs(base, S_sing[keep, , drop = FALSE],
                                   sing_labels[keep], S_dbl, dbl_pair,
                                   lib_sing[keep], lambda = 0.5, R = 300L,
                                   seed = 3L)
  for (pl in grep("D", names(base$pair_models), value = TRUE)) {
    expect_identical(out$models$pair_models[[pl]], base$pair_models[[pl]])
  }
  untouched <- setdiff(names(base$pair_models), grep("D", names(base$pair_models), value = TRUE))
  expect_false(isTRUE(all.equal(out$models$pair_models[[untouched[1]]],
                                base$pair_models[[untouched[1]]])))
})

test_that("shrunk models remain valid multivariate normals", {
  out <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels, S_dbl, dbl_pair,
                                   lib_sing, lambda = 0.5, R = 300L, seed = 3L)
  for (pl in names(out$models$pair_models)) {
    S <- recover_Sigma(out$models$pair_models[[pl]])
    expect_true(isSymmetric(S, tol = 1e-8))
    expect_true(all(eigen(S, symmetric = TRUE, only.values = TRUE)$values > 0))
    expect_true(is.finite(out$models$pair_models[[pl]]$log_det))
  }
  ll <- Obligato:::all_pair_ll(S_dbl[1:20, , drop = FALSE],
                               out$models$pair_models)
  expect_true(all(is.finite(ll)))
})

test_that("too little data to cross-validate falls back to no shrinkage", {
  # Every pair below min_pair, so no fold can be scored. The procedure must
  # decline to shrink rather than fail or choose arbitrarily.
  few <- unlist(lapply(pairs_all, function(p) head(which(dbl_pair ==
    Obligato:::pair_label(p[1], p[2])), 4)))
  out <- Obligato:::apply_eb_pairs(base, S_sing, sing_labels,
                                   S_dbl[few, , drop = FALSE], dbl_pair[few],
                                   lib_sing, lambda = "cv", grid = c(0, 0.5),
                                   R = 200L, seed = 3L)
  expect_identical(out$lambda, 0)
  expect_identical(out$models$pair_models, base$pair_models)
})
