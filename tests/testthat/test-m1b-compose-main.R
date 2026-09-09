################################################################################
# M1b: the high-level compose() orchestrator. Its faithful end-to-end
# reproduction of the King/batch2 numbers is the M2b integration gate (needs
# PhiSpace + real data). Here we test the ORCHESTRATION with an injected
# deterministic scorer: the detection module branches, the returned structure,
# determinism, and that compose() invokes compose_pairs() correctly.
################################################################################

# Deterministic synthetic scorer: cells x types scores in (-0.9, 0.9), depending
# only on position, so it never touches the global RNG stream.
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
  reference <- list(celltype = rep(types, length.out = 300))
  hc <- rep(FALSE, 200); hc[1:120] <- TRUE
  hc_labels <- rep(types, each = 20)   # 120 hc singlets, 20 per type
  list(query = query, reference = reference, hc = hc, hc_labels = hc_labels, types = types)
}

common_args <- function(s, ...) {
  list(query = s$query, reference = s$reference, phenotypes = "celltype",
       hc_singlets = s$hc, hc_labels = s$hc_labels, n_per_pair = 30,
       score_fn = fake_scorer, ...)
}

test_that("compose builtin detector returns composition + score + flag", {
  skip_if_not_installed("xgboost")
  s <- setup_query()
  res <- do.call(compose, common_args(s, detector = "builtin"))
  expect_s3_class(res, "obligato_composition")
  expect_equal(nrow(res$composition), 200)
  expect_true(all(c("cell", "c1", "c2") %in% names(res$composition)))
  expect_length(res$detection_score, 200)
  expect_length(res$flag, 200)
  expect_true(is.logical(res$flag))
  expect_true(is.numeric(res$info$threshold))
  expect_false(is.null(res$detection_model))
  expect_identical(res$info$detection_features, "PhiSpace scores plus library size")
})

test_that("detector = 'none' composes without flagging", {
  s <- setup_query()
  res <- do.call(compose, common_args(s, detector = "none"))
  expect_null(res$detection_score)
  expect_null(res$flag)
  expect_equal(nrow(res$composition), 200)
})

test_that("bring-your-own-detector uses the supplied flags", {
  s <- setup_query()
  my_flags <- rep(FALSE, 200); my_flags[c(5, 50, 150)] <- TRUE
  res <- do.call(compose, common_args(s, detector = my_flags))
  expect_identical(res$flag, my_flags)
  expect_null(res$detection_score)
  expect_identical(res$info$detector, "byo")

  # Numeric scores require an explicit threshold and are retained in the result.
  sc <- runif(200)
  expect_error(do.call(compose, common_args(s, detector = sc)), "detector_threshold")
  res2 <- do.call(compose, common_args(
    s, detector = sc, detector_threshold = 0.5))
  expect_identical(res2$flag, sc > 0.5)
})

test_that("compose is deterministic and invokes compose_pairs correctly", {
  skip_if_not_installed("xgboost")
  s <- setup_query()
  r1 <- do.call(compose, common_args(s, detector = "builtin"))
  r2 <- do.call(compose, common_args(s, detector = "builtin"))
  expect_identical(r1$composition, r2$composition)
  expect_identical(r1$flag, r2$flag)

  # composition equals compose_pairs() run on the same scores + fitted models
  S_all <- fake_scorer(s$query, s$reference, "celltype")
  manual <- compose_pairs(S_all, r1$models$pair_models, candidates = "all5", output = "top1")
  expect_identical(r1$composition$c1, manual$c1)
  expect_identical(r1$composition$c2, manual$c2)
})

test_that("an invalid detector argument errors clearly", {
  s <- setup_query()
  expect_error(do.call(compose, common_args(s, detector = "magic")),
               "must be")
})
