################################################################################
# Same-type composition (D2) and the frequency pair prior, adopted 2026-09-29.
# Evidence and design: project paper_plan.md §4.1g. The contract tested here:
#   * same-type training doublets sum two DISTINCT singlets of one type;
#   * compose_pairs() adds exactly one same-type candidate, {top1, top1};
#   * the prior changes ranking only, never the reported log-likelihood;
#   * with both options off, composition equals the heterotypic-only ranker;
#   * detection is identical with same_type on or off.
################################################################################

fake_scorer <- function(counts, reference, phenotypes, label = "") {
  types <- sort(unique(as.character(reference[[phenotypes]])))
  m <- ncol(counts); K <- length(types)
  M <- outer(seq_len(m), seq_len(K), function(i, j) 0.9 * sin(0.3 * i + 0.7 * j + 0.05 * i * j))
  colnames(M) <- types
  M
}
setup_query <- function() {
  types <- c("A", "B", "C", "D", "E", "F")
  set.seed(11)
  query <- matrix(rpois(40 * 200, 3), nrow = 40,
                  dimnames = list(paste0("g", 1:40), paste0("cell", 1:200)))
  list(query = query, reference = list(celltype = rep(types, length.out = 300)),
       hc = c(rep(TRUE, 120), rep(FALSE, 80)), hc_labels = rep(types, each = 20))
}
args_for <- function(s, ...) list(query = s$query, reference = s$reference,
  phenotypes = "celltype", hc_singlets = s$hc, hc_labels = s$hc_labels,
  n_per_pair = 30, score_fn = fake_scorer, ...)

test_that("same-type doublets sum two distinct cells of one type, reproducibly", {
  cnt <- matrix(seq_len(5 * 7), 5, 7)
  hc_idx <- 1:7; hc_types <- c("A", "A", "A", "B", "B", "C", "D")
  before <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  sim <- Obligato:::simulate_same_type_doublets(cnt, hc_idx, hc_types, n_per_pair = 40, seed = 3)
  after <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  expect_identical(before, after)
  expect_setequal(unique(sim$pair_label), c("A + A", "B + B"))   # C and D have one cell
  expect_true(all(sim$idx_A != sim$idx_B))
  expect_true(all(hc_types[sim$idx_A] == hc_types[sim$idx_B]))
  expect_equal(unname(sim$counts), unname(cnt[, sim$idx_A] + cnt[, sim$idx_B]))
  sim2 <- Obligato:::simulate_same_type_doublets(cnt, hc_idx, hc_types, n_per_pair = 40, seed = 3)
  expect_identical(sim, sim2)
})

test_that("pair log prior is the two-draw probability", {
  p <- c(A = 0.5, B = 0.3, C = 0.2)
  lp <- Obligato:::pair_log_prior(c("A + B", "A + A", "B + C"), p)
  expect_equal(unname(exp(lp)), c(2 * 0.5 * 0.3, 0.25, 2 * 0.3 * 0.2))
  expect_true(is.na(Obligato:::pair_log_prior("A + Z", p)))
})

test_that("compose_pairs adds exactly {top1,top1} and matches a brute-force oracle", {
  set.seed(5)
  types <- LETTERS[1:6]
  S <- matrix(runif(60 * 6, -0.9, 0.9), 60, 6, dimnames = list(NULL, types))
  mk <- function() { mu <- rnorm(6); list(mu = mu, Sigma_inv = diag(6), log_det = 0) }
  het <- setNames(lapply(combn(types, 2, simplify = FALSE), function(x) mk()),
                  vapply(combn(types, 2, simplify = FALSE), function(x) paste(x, collapse = " + "), ""))
  same <- setNames(lapply(types, function(x) mk()), paste(types, "+", types))
  lp <- Obligato:::pair_log_prior(c(names(het), names(same)), setNames(c(.3, .2, .2, .1, .1, .1), types))
  for (prior in list(NULL, lp)) {
    got <- compose_pairs(S, het, same_models = same, log_prior = prior)
    Y <- Obligato:::to_logit(S)
    for (i in seq_len(nrow(S))) {
      top <- types[order(S[i, ], decreasing = TRUE)[1:5]]
      cand <- c(combn(top, 2, FUN = function(x) Obligato:::pair_label(x[1], x[2])),
                Obligato:::pair_label(top[1], top[1]))
      sc <- vapply(cand, function(pl) Obligato:::mvn_ll_rows(Y[i, , drop = FALSE],
                   c(het, same)[[pl]]), 0)
      if (!is.null(prior)) sc <- sc + prior[cand]
      expect_identical(Obligato:::pair_label(got$c1[i], got$c2[i]), cand[which.max(sc)])
    }
  }
  # the same-type candidate is the TOP type only
  diag_calls <- got$c1 == got$c2
  top1 <- types[max.col(S, ties.method = "first")]
  expect_true(all(got$c1[diag_calls] == top1[diag_calls]))
})

test_that("topk keeps the raw log-likelihood and reports the prior separately", {
  set.seed(6)
  types <- LETTERS[1:4]
  S <- matrix(runif(20 * 4, -0.9, 0.9), 20, 4, dimnames = list(NULL, types))
  mk <- function() list(mu = rnorm(4), Sigma_inv = diag(4), log_det = 0)
  het <- setNames(lapply(1:6, function(i) mk()),
                  vapply(combn(types, 2, simplify = FALSE), function(x) paste(x, collapse = " + "), ""))
  same <- setNames(lapply(types, function(x) mk()), paste(types, "+", types))
  lp <- Obligato:::pair_log_prior(c(names(het), names(same)), setNames(rep(.25, 4), types))
  tk <- compose_pairs(S, het, output = "topk", k = 3, same_models = same, log_prior = lp)
  expect_true("log_prior" %in% names(tk))
  Y <- Obligato:::to_logit(S)
  raw <- vapply(seq_len(nrow(tk)), function(r) Obligato:::mvn_ll_rows(
    Y[tk$cell[r], , drop = FALSE], c(het, same)[[Obligato:::pair_label(tk$c1[r], tk$c2[r])]]), 0)
  expect_equal(tk$ll, raw)
  by_cell <- split(tk$ll + tk$log_prior, tk$cell)
  expect_true(all(vapply(by_cell, function(v) all(diff(v) <= 1e-12), TRUE)))
  expect_error(compose_pairs(S, het, same_models = same, log_prior = lp[1:3]), "covering every model")
  expect_error(compose_pairs(S, het, same_models = het), "same-type pairs only")
})

test_that("with both options off, composition equals the heterotypic-only ranker", {
  s <- setup_query()
  res <- do.call(compose, args_for(s, detector = "none", same_type = FALSE, pair_prior = "none"))
  S_all <- fake_scorer(s$query, s$reference, "celltype")
  ref <- compose_pairs(S_all, res$models$pair_models, candidates = "all5", output = "top1")
  expect_identical(res$composition, ref)
  expect_null(res$models$same_models)
  expect_false(res$info$same_type)
})

test_that("detection is identical with same_type on or off", {
  skip_if_not_installed("xgboost")
  s <- setup_query()
  on <- do.call(compose, args_for(s, detector = "builtin"))
  off <- do.call(compose, args_for(s, detector = "builtin", same_type = FALSE, pair_prior = "none"))
  expect_identical(on$detection_score, off$detection_score)
  expect_identical(on$flag, off$flag)
  expect_identical(on$detection_features, off$detection_features)
  expect_identical(on$models$pair_models, off$models$pair_models)
  expect_true(on$info$same_type)
  expect_identical(on$info$pair_prior, "frequency")
  expect_gt(on$info$n_same_models, 0)
  expect_equal(sum(on$info$type_freq), 1)
})
