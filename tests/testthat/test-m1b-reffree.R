################################################################################
# M2c: reference-free mode compose_reffree(). Tested with injected deterministic
# clustering / normalisation / scoring so no Seurat or PhiSpace is required — the
# loop's orchestration, the ecdf flag, the cluster-pair composition, and the
# returned structure are exercised end-to-end.
################################################################################
library(SingleCellExperiment)

# 4 fixed clusters by position -> guarantees >= sing_min cells/cluster and >= min_pair
# doublets/pair so every model fits.
fake_cluster <- function(counts, resolution = 1.0, npcs = 30L, nfeatures = 2000L)
  paste0("c", (seq_len(ncol(counts)) %% 4L))

# Add a logcounts assay without scran.
fake_norm <- function(sce) { logcounts(sce) <- log1p(as.matrix(counts(sce))); sce }

# Return a cells x clusters score matrix in (-0.9, 0.9), colnames = cluster ids,
# depending only on position so it never touches the RNG stream.
fake_score <- function(counts, ref_sce, phenotypes, label = "") {
  types <- sort(unique(as.character(ref_sce[[phenotypes]])))
  m <- ncol(counts); K <- length(types)
  M <- outer(seq_len(m), seq_len(K), function(i, j) 0.9 * sin(0.25 * i + 0.6 * j))
  colnames(M) <- types
  M
}

setup_query <- function(n = 240L, g = 30L) {
  set.seed(3)
  q <- matrix(rpois(g * n, 3), nrow = g,
              dimnames = list(paste0("g", 1:g), paste0("cell", 1:n)))
  q
}

run <- function(...) compose_reffree(
  query = setup_query(), target_rate = 0.1, n_rounds = 1L, min_pair = 6L, sing_min = 5L,
  cluster_fn = fake_cluster, score_fn = fake_score, norm_fn = fake_norm, seed = 1L, ...)

test_that("compose_reffree returns a well-formed obligato_reffree object", {
  skip_if_not_installed("xgboost")
  res <- run()
  expect_s3_class(res, "obligato_reffree")
  N <- ncol(setup_query())
  expect_equal(nrow(res$composition), N)
  expect_true(all(c("cell","top_pair","top_loglik") %in% names(res$composition)))
  expect_length(res$detection_score, N)
  expect_length(res$flag, N); expect_true(is.logical(res$flag))
  expect_length(res$rounds, 2L)                       # rounds 0 and 1
  expect_length(res$clusters, N)
  expect_identical(res$info$detection_features, "PhiSpace scores")
  expect_true(length(res$pair_models) > 0L)
  expect_true(all(res$composition$top_pair %in% names(res$pair_models)))
})

test_that("flag uses the per-round ecdf threshold, not 0.5", {
  skip_if_not_installed("xgboost")
  res <- run()
  expect_identical(res$flag, res$detection_score > res$info$threshold)
  # the info threshold is the final round's ecdf crossover
  expect_equal(res$info$threshold, res$rounds[[length(res$rounds)]]$threshold)
})

test_that("compose_reffree is deterministic given seed", {
  skip_if_not_installed("xgboost")
  r1 <- run(); r2 <- run()
  expect_identical(r1$detection_score, r2$detection_score)
  expect_identical(r1$composition$top_pair, r2$composition$top_pair)
  expect_identical(r1$info$threshold, r2$info$threshold)
})

test_that("hc_singlets seeds round 0 (its complement is the initial doublet flag)", {
  skip_if_not_installed("xgboost")
  N <- ncol(setup_query())
  hc <- rep(TRUE, N); hc[1:20] <- FALSE          # 20 seeded doublets
  res <- compose_reffree(query = setup_query(), hc_singlets = hc, target_rate = 0.1,
                         n_rounds = 0L, min_pair = 6L, sing_min = 5L,
                         cluster_fn = fake_cluster, score_fn = fake_score,
                         norm_fn = fake_norm, seed = 1L)
  # round 0 clustered the 220 seeded singlets, not all 240
  expect_equal(res$rounds[[1]]$n_clean, N - 20L)
})

test_that("target_rate is validated", {
  expect_error(compose_reffree(setup_query(), target_rate = 0,  cluster_fn = fake_cluster,
                               score_fn = fake_score, norm_fn = fake_norm), "target_rate")
  expect_error(compose_reffree(setup_query(), target_rate = 1.5, cluster_fn = fake_cluster,
                               score_fn = fake_score, norm_fn = fake_norm), "target_rate")
})

test_that("simulate_cluster_doublets floors per-pair count at min_pair", {
  set.seed(1)
  cnt <- matrix(rpois(20 * 200, 3), nrow = 20)
  cl  <- paste0("c", seq_len(200) %% 5L)          # 5 clusters -> 10 pairs
  sim <- simulate_cluster_doublets(cnt, seq_len(200), cl, target_rate = 0.02, min_pair = 6L, seed = 1L)
  # low target rate -> ceiling < min_pair, so each realised pair has exactly min_pair sims
  expect_true(all(table(sim$pair_label) >= 6L))
  expect_equal(ncol(sim$counts), length(sim$pair_label))
})
