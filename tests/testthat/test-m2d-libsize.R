# M2d — library-size detection features and multi-level scoring (2026-09-09).
# Exactness and regression coverage for the production M2d orchestration.

test_that("libsize_features computes log_total, n_genes, typerel exactly", {
  set.seed(1)
  cts <- matrix(rpois(30 * 8, 4), 30, 8); cts[1:6, 1] <- 0   # gene x cell
  S <- matrix(runif(8 * 3), 8, 3); colnames(S) <- c("A", "B", "C")
  top1 <- colnames(S)[max.col(S, ties.method = "first")]
  tm <- c(A = log1p(50), B = log1p(60), C = log1p(40))

  lf <- Obligato:::libsize_features(cts, top1, tm)
  expect_equal(lf$log_total, log1p(Matrix::colSums(cts)))
  expect_equal(lf$n_genes,   as.numeric(Matrix::colSums(cts > 0)))
  expect_equal(lf$typerel,   unname(log1p(Matrix::colSums(cts)) - tm[top1]))
  expect_setequal(names(lf), c("log_total", "n_genes", "typerel"))
})

test_that("libsize_features omits typerel when no type_median (raw-only, A1)", {
  cts <- matrix(rpois(30 * 4, 4), 30, 4)
  lf <- Obligato:::libsize_features(cts, rep("A", 4), type_median = NULL)
  expect_setequal(names(lf), c("log_total", "n_genes"))
})

test_that("typerel is 0 when the top-1 type has no median (unknown type)", {
  cts <- matrix(rpois(30 * 4, 4), 30, 4)
  lf <- Obligato:::libsize_features(cts, rep("Z", 4), type_median = c(A = 3))
  expect_true(all(lf$typerel == 0))
})

test_that("type-relative medians use score-derived top-1 types", {
  cts <- matrix(c(4, 5, 9, 10, 49, 50, 54, 55), nrow = 2)
  S <- rbind(c(0.9, 0.1), c(0.1, 0.9), c(0.8, 0.2), c(0.2, 0.8))
  colnames(S) <- c("A", "B")

  got <- Obligato:::libsize_type_medians(cts, S)
  expect_equal(unname(got["A"]), stats::median(log1p(c(9, 99))))
  expect_equal(unname(got["B"]), stats::median(log1p(c(19, 109))))
  expect_error(Obligato:::libsize_type_medians(cts[, 1:3], S), "must match")
})

test_that("compute_features suffixes feature names for concatenated levels", {
  set.seed(2)
  S <- matrix(runif(10 * 4), 10, 4); colnames(S) <- LETTERS[1:4]
  Y <- Obligato:::to_logit(S)
  sm <- list(A = Obligato:::fit_ln_shrink(Y[1:5, , drop = FALSE]))
  pm <- list(`A + B` = Obligato:::fit_ln_shrink(Y[1:6, , drop = FALSE]))

  f0 <- Obligato:::compute_features(S, sm, pm)                 # unsuffixed (single level)
  f2 <- Obligato:::compute_features(S, sm, pm, suffix = "_l2") # suffixed
  expect_equal(ncol(f0), 7)
  expect_equal(names(f2), paste0(names(f0), "_l2"))
  expect_equal(unname(as.matrix(f2)), unname(as.matrix(f0)))   # values unchanged
})

# --- multi-level detection (M2d option ii): concatenate L1 + L2 features --------
.ml_fake_scorer <- function(counts, reference, phenotypes, label = "") {
  types <- sort(unique(as.character(reference[[phenotypes]])))
  m <- ncol(counts); K <- length(types)
  M <- outer(seq_len(m), seq_len(K), function(i, j) 0.9 * sin(0.3 * i + 0.7 * j + 0.05 * i * j))
  colnames(M) <- types; M
}
.ml_setup <- function() {
  types <- c("A", "B", "C", "D", "E", "F")
  set.seed(11)
  query <- matrix(rpois(40 * 200, 3), nrow = 40,
                  dimnames = list(paste0("g", 1:40), paste0("cell", 1:200)))
  ref1  <- list(celltype  = rep(types, length.out = 300))
  ref2  <- list(celltype2 = rep(c("X", "Y", "Z"), length.out = 300))  # coarser 2nd level
  hc <- rep(FALSE, 200); hc[1:120] <- TRUE
  hc_labels <- rep(types, each = 20)
  list(query = query, ref1 = ref1, ref2 = ref2, hc = hc, hc_labels = hc_labels, types = types)
}

test_that("compose runs multi-level detection and composition uses the primary level", {
  skip_if_not_installed("xgboost")
  s <- .ml_setup()
  res <- compose(query = s$query, reference = list(s$ref1, s$ref2),
                 phenotypes = c("celltype", "celltype2"),
                 hc_singlets = s$hc, hc_labels = s$hc_labels, n_per_pair = 30,
                 detector = "builtin", score_fn = .ml_fake_scorer)
  expect_s3_class(res, "obligato_composition")
  expect_length(res$detection_score, 200)
  # composition is driven by the PRIMARY level only (types A..F), not the 2nd level
  expect_true(all(res$composition$c1 %in% s$types))
})

test_that("adding a 2nd annotation level changes the detector (features concatenated)", {
  skip_if_not_installed("xgboost")
  s <- .ml_setup()
  one <- compose(query = s$query, reference = s$ref1, phenotypes = "celltype",
                 hc_singlets = s$hc, hc_labels = s$hc_labels, n_per_pair = 30,
                 detector = "builtin", score_fn = .ml_fake_scorer)
  two <- compose(query = s$query, reference = list(s$ref1, s$ref2),
                 phenotypes = c("celltype", "celltype2"),
                 hc_singlets = s$hc, hc_labels = s$hc_labels, n_per_pair = 30,
                 detector = "builtin", score_fn = .ml_fake_scorer)
  # same composition (primary level unchanged), different detection scores
  expect_equal(one$composition$c1, two$composition$c1)
  expect_false(isTRUE(all.equal(one$detection_score, two$detection_score)))
})

test_that("multi-level conformal composition stays in the primary score space", {
  s <- .ml_setup()
  res <- compose(query = s$query, reference = list(s$ref1, s$ref2),
                 phenotypes = c("celltype", "celltype2"),
                 hc_singlets = s$hc, hc_labels = s$hc_labels, n_per_pair = 30,
                 detector = "none", output = "conformal", alpha = 0.10,
                 score_fn = .ml_fake_scorer)

  expect_s3_class(res, "obligato_composition")
  expect_length(res$conformal$sizes, ncol(s$query))
  expect_true(all(res$composition$c1 %in% s$types))
  expect_true(all(res$composition$c2 %in% s$types))
})
