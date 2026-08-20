################################################################################
# M1b GATE: fit_doublet_models() must reproduce Section 3 of
# run_compose_kingetal.R:227-247 exactly (identical block in
# run_methodC_compose.R's process_batch fitting loop).
#
# The block is inline, not a function, so its arithmetic is transcribed here as an
# independent oracle and the extracted function is required identical to it.
################################################################################

# Independent transcription of the script's Section 3 loops.
oracle_fit <- function(S_train_sing, hc_types, S_train_dbl, dbl_pl, types,
                       min_singlet = 15, min_pair = 30) {
  Y_ts <- to_logit(S_train_sing)
  sing_models <- list()
  for (tt in types) {
    idx <- which(hc_types == tt)
    if (length(idx) >= min_singlet) sing_models[[tt]] <- fit_ln_shrink(Y_ts[idx, , drop = FALSE])
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
    if (length(idx) >= 30) pair_models[[pl]] <- fit_ln_shrink(Y_td[idx, , drop = FALSE])
  }
  list(sing_models = sing_models, sing_marginal_sd = sing_marginal_sd, pair_models = pair_models)
}

# Build a deterministic training pool: `types`, `n_sing_per` singlets each,
# all heterotypic pairs, `n_dbl_per` doublets each. Scores live in (-1, 1).
make_pool <- function(types, n_sing_per, n_dbl_per, seed = 1) {
  set.seed(seed)
  K <- length(types)
  hc_types <- rep(types, times = n_sing_per)
  S_sing <- matrix(runif(length(hc_types) * K, -0.98, 0.98), ncol = K,
                   dimnames = list(NULL, types))
  prs <- combn(types, 2, simplify = FALSE)
  dbl_pl <- unlist(lapply(prs, function(p) rep(pair_label(p[1], p[2]), n_dbl_per)))
  S_dbl <- matrix(runif(length(dbl_pl) * K, -0.98, 0.98), ncol = K,
                  dimnames = list(NULL, types))
  list(S_sing = S_sing, hc_types = hc_types, S_dbl = S_dbl, dbl_pl = dbl_pl)
}

test_that("fit_doublet_models is identical to the script's Section 3", {
  types <- c("A", "B", "C", "D")
  p <- make_pool(types, n_sing_per = c(A = 25, B = 25, C = 25, D = 25), n_dbl_per = 40)
  got <- fit_doublet_models(p$S_sing, p$hc_types, p$S_dbl, p$dbl_pl, types = types)
  exp <- oracle_fit(p$S_sing, p$hc_types, p$S_dbl, p$dbl_pl, types = types)
  expect_equal(got, exp, tolerance = 0)
})

test_that("min_singlet and min_pair cutoffs drop under-supported models exactly as the script", {
  types <- c("A", "B", "C", "D")
  # A has only 10 singlets (< 15) -> dropped; others have 25.
  p <- make_pool(types, n_sing_per = c(A = 10, B = 25, C = 25, D = 25), n_dbl_per = 40)
  # thin one pair below 30 doublets
  keep <- !(p$dbl_pl == pair_label("A", "B") & ave(seq_along(p$dbl_pl), p$dbl_pl, FUN = seq_along) > 20)
  S_dbl <- p$S_dbl[keep, , drop = FALSE]; dbl_pl <- p$dbl_pl[keep]

  got <- fit_doublet_models(p$S_sing, p$hc_types, S_dbl, dbl_pl, types = types)
  exp <- oracle_fit(p$S_sing, p$hc_types, S_dbl, dbl_pl, types = types)
  expect_equal(got, exp, tolerance = 0)

  expect_false("A" %in% names(got$sing_models))          # 10 < 15
  expect_true(all(c("B", "C", "D") %in% names(got$sing_models)))
  expect_false(pair_label("A", "B") %in% names(got$pair_models))  # 20 < 30
  expect_identical(names(got$sing_models), names(got$sing_marginal_sd))
})

test_that("type iteration order changes only list order, never a fitted model", {
  types <- c("A", "B", "C", "D")
  p <- make_pool(types, n_sing_per = c(A = 25, B = 25, C = 25, D = 25), n_dbl_per = 40)
  fwd <- fit_doublet_models(p$S_sing, p$hc_types, p$S_dbl, p$dbl_pl, types = types)
  rev <- fit_doublet_models(p$S_sing, p$hc_types, p$S_dbl, p$dbl_pl, types = rev(types))
  expect_setequal(names(fwd$sing_models), names(rev$sing_models))
  for (tt in names(fwd$sing_models)) expect_equal(fwd$sing_models[[tt]], rev$sing_models[[tt]], tolerance = 0)
})

test_that("default types equals sorted unique singlet types", {
  types <- c("D", "A", "C", "B")
  p <- make_pool(types, n_sing_per = c(D = 25, A = 25, C = 25, B = 25), n_dbl_per = 40)
  got_default <- fit_doublet_models(p$S_sing, p$hc_types, p$S_dbl, p$dbl_pl)
  got_sorted  <- fit_doublet_models(p$S_sing, p$hc_types, p$S_dbl, p$dbl_pl, types = sort(unique(p$hc_types)))
  expect_equal(got_default, got_sorted, tolerance = 0)
})
