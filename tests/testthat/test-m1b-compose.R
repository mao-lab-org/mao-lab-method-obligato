################################################################################
# M1b GATE: compose() must reproduce Section 6 of run_compose_kingetal.R:342-397
# (the DblLik branch), identical to run_methodC_compose.R's composition block.
#
# The block is inline, so its loop is transcribed here as an oracle. compose(...,
# output = "top1") must return the same c1/c2 as the script for every candidate
# set. output = "topk" and the bring-your-own-detector `flagged` path are new
# surface tested against the same ranking.
################################################################################

# Exact transcription of the script's Section 6 DblLik loop for one candidate set.
oracle_top1 <- function(S, pair_models, cs) {
  n <- nrow(S); dims <- colnames(S); Y <- to_logit(S)
  topK <- t(apply(S, 1, function(r) dims[order(r, decreasing = TRUE)[1:5]]))
  c1 <- character(n); c2 <- character(n)
  for (i in seq_len(n)) {
    cand_df <- make_candidates(topK[i, ])[[cs]]
    best_c1 <- cand_df$c1[1]; best_c2 <- cand_df$c2[1]; best_ll <- -Inf
    for (j in seq_len(nrow(cand_df))) {
      pl <- pair_label(cand_df$c1[j], cand_df$c2[j])
      if (pl %in% names(pair_models)) {
        ll <- mvn_ll_rows(matrix(Y[i, ], nrow = 1), pair_models[[pl]])
        if (ll > best_ll) { best_ll <- ll; best_c1 <- cand_df$c1[j]; best_c2 <- cand_df$c2[j] }
      }
    }
    c1[i] <- best_c1; c2[i] <- best_c2
  }
  data.frame(c1 = c1, c2 = c2, stringsAsFactors = FALSE)
}

# Build fitted pair models and a fresh query score matrix.
setup_scenario <- function(seed = 7) {
  types <- c("A", "B", "C", "D", "E", "F")
  K <- length(types)
  set.seed(seed)
  hc_types <- rep(types, each = 30)
  S_sing <- matrix(runif(length(hc_types) * K, -0.98, 0.98), ncol = K, dimnames = list(NULL, types))
  prs <- combn(types, 2, simplify = FALSE)
  dbl_pl <- unlist(lapply(prs, function(p) rep(pair_label(p[1], p[2]), 40)))
  S_dbl  <- matrix(runif(length(dbl_pl) * K, -0.98, 0.98), ncol = K, dimnames = list(NULL, types))
  mods <- fit_doublet_models(S_sing, hc_types, S_dbl, dbl_pl, types = types)
  S_test <- matrix(runif(50 * K, -0.98, 0.98), ncol = K, dimnames = list(NULL, types))
  list(pair_models = mods$pair_models, S = S_test, types = types)
}

test_that("compose top1 reproduces Section 6 for all three candidate sets", {
  sc <- setup_scenario()
  for (cs in c("all5", "flex2", "fix1")) {
    got <- compose_pairs(sc$S, sc$pair_models, candidates = cs, output = "top1")
    exp <- oracle_top1(sc$S, sc$pair_models, cs)
    expect_identical(got$c1, exp$c1, info = cs)
    expect_identical(got$c2, exp$c2, info = cs)
  }
  # default candidate set is the locked all5
  expect_identical(compose_pairs(sc$S, sc$pair_models)$c1, oracle_top1(sc$S, sc$pair_models, "all5")$c1)
})

test_that("compose top1 still matches the oracle when a pair model is absent", {
  sc <- setup_scenario()
  drop <- names(sc$pair_models)[1]
  pm <- sc$pair_models[setdiff(names(sc$pair_models), drop)]
  got <- compose_pairs(sc$S, pm, candidates = "all5", output = "top1")
  exp <- oracle_top1(sc$S, pm, "all5")
  expect_identical(got$c1, exp$c1)
  expect_identical(got$c2, exp$c2)
})

test_that("topk rank 1 equals top1, and ranks are ordered by likelihood", {
  sc <- setup_scenario()
  top1 <- compose_pairs(sc$S, sc$pair_models, output = "top1")
  topk <- compose_pairs(sc$S, sc$pair_models, output = "topk", k = 3)
  r1 <- topk[topk$rank == 1, ]
  expect_identical(r1$c1, top1$c1)
  expect_identical(r1$c2, top1$c2)
  # within each cell, ll is non-increasing with rank
  for (cc in unique(topk$cell)) {
    lls <- topk$ll[topk$cell == cc]
    expect_false(is.unsorted(rev(lls)))
  }
  expect_lte(max(table(topk$cell)), 3)
})

test_that("flagged composes only the flagged cells, identically", {
  sc <- setup_scenario()
  flag <- rep(FALSE, nrow(sc$S)); flag[c(3, 10, 25, 41)] <- TRUE
  got <- compose_pairs(sc$S, sc$pair_models, flagged = flag, output = "top1")
  full <- oracle_top1(sc$S, sc$pair_models, "all5")
  expect_identical(got$cell, which(flag))
  expect_identical(got$c1, full$c1[which(flag)])
  expect_identical(got$c2, full$c2[which(flag)])
})

test_that("conformal output errors clearly (deferred to its own step)", {
  sc <- setup_scenario()
  expect_error(compose_pairs(sc$S, sc$pair_models, output = "conformal"), "not yet implemented")
})
