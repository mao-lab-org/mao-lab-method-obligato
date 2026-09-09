.uf_scorer <- function(counts, reference, phenotypes, label = "") {
  types <- sort(unique(as.character(reference[[phenotypes]])))
  M <- outer(seq_len(ncol(counts)), seq_along(types),
             function(i, j) 0.9 * sin(0.3 * i + 0.7 * j + 0.05 * i * j))
  colnames(M) <- types
  M
}

.uf_setup <- function() {
  types <- LETTERS[1:6]
  set.seed(25)
  query <- matrix(rpois(40 * 200, 3), 40, 200,
                  dimnames = list(paste0("g", 1:40), paste0("c", 1:200)))
  reference <- list(celltype = rep(types, 50))
  hc <- rep(FALSE, 200); hc[1:120] <- TRUE
  labels <- rep(types, each = 20)
  list(query = query, reference = reference, hc = hc, labels = labels)
}

.uf_args <- function(s, ...) {
  list(query = s$query, reference = s$reference, phenotypes = "celltype",
       hc_singlets = s$hc, hc_labels = s$labels, n_per_pair = 30,
       score_fn = .uf_scorer, ...)
}

test_that("full-length hc_labels are subset to high-confidence singlets", {
  types <- LETTERS[1:6]
  set.seed(21)
  query <- matrix(rpois(30 * 180, 3), 30, 180,
                  dimnames = list(paste0("g", 1:30), paste0("c", 1:180)))
  reference <- list(celltype = rep(types, 40))
  hc <- rep(FALSE, 180); hc[1:120] <- TRUE
  labels <- rep(types, length.out = 180)
  scorer <- function(counts, reference, phenotypes, label = "") {
    ty <- sort(unique(as.character(reference[[phenotypes]])))
    z <- matrix(runif(ncol(counts) * length(ty), -0.8, 0.8),
                nrow = ncol(counts), dimnames = list(NULL, ty))
    z
  }

  res <- compose(query, reference, "celltype", hc, hc_labels = labels,
                 detector = "none", n_per_pair = 30, score_fn = scorer)
  expect_identical(res$info$n_hc_singlets, 120L)
  expect_false(anyNA(res$composition$c1))
})

test_that("compose obeys seed without changing the caller RNG state", {
  types <- LETTERS[1:5]
  set.seed(22)
  query <- matrix(rpois(20 * 100, 3), 20, 100,
                  dimnames = list(paste0("g", 1:20), paste0("c", 1:100)))
  reference <- list(celltype = rep(types, 30))
  hc <- rep(TRUE, 100)
  labels <- rep(types, each = 20)
  scorer <- function(counts, reference, phenotypes, label = "") {
    ty <- sort(unique(as.character(reference[[phenotypes]])))
    matrix(runif(ncol(counts) * length(ty), -0.8, 0.8),
           nrow = ncol(counts), dimnames = list(NULL, ty))
  }

  set.seed(987)
  before <- .Random.seed
  r1 <- compose(query, reference, "celltype", hc, labels,
                detector = "none", n_per_pair = 30, seed = 4, score_fn = scorer)
  expect_identical(.Random.seed, before)
  r2 <- compose(query, reference, "celltype", hc, labels,
                detector = "none", n_per_pair = 30, seed = 4, score_fn = scorer)
  expect_identical(r1$composition, r2$composition)
})

test_that("composition supports two-to-four type references and stable empty top-k", {
  set.seed(23)
  S <- matrix(runif(40, -0.8, 0.8), 10, 4,
              dimnames = list(NULL, LETTERS[1:4]))
  Y <- Obligato:::to_logit(S)
  prs <- combn(colnames(S), 2, simplify = FALSE)
  pm <- lapply(prs, function(x) Obligato:::fit_ln_shrink(rbind(Y, Y + 0.01)))
  names(pm) <- vapply(prs, function(x) Obligato:::pair_label(x[1], x[2]), character(1))

  got <- Obligato:::compose_pairs(S, pm)
  expect_equal(nrow(got), nrow(S))
  expect_false(anyNA(got$c1))
  empty <- Obligato:::compose_pairs(
    S, pm, flagged = rep(FALSE, nrow(S)), output = "topk")
  expect_s3_class(empty, "data.frame")
  expect_identical(names(empty), c("cell", "rank", "c1", "c2", "ll"))
  expect_equal(nrow(empty), 0L)
})

test_that("conformal calibration uses the finite-sample order statistic", {
  nonconformity <- 1:10
  log_h <- matrix(-nonconformity, ncol = 1,
                  dimnames = list(NULL, "A + B"))
  q <- Obligato:::conformal_calibrate(
    log_h, rep("A + B", 10), alpha = 0.20)[[1]]
  expect_equal(q, 9)
})

test_that("AUPRC is invariant to row order when scores are tied", {
  truth <- c(1, 0, 1, 0, 1, 0)
  scores <- c(0.8, 0.8, 0.5, 0.5, 0.2, 0.2)
  set.seed(24)
  perm <- sample(seq_along(truth))
  a <- compute_auprc(truth, scores)
  b <- compute_auprc(truth[perm], scores[perm])
  expect_equal(a$auprc, b$auprc)
  expect_equal(a$curve, b$curve)
})

test_that("results carry user-facing provenance and fitted detector state", {
  skip_if_not_installed("xgboost")
  s <- .uf_setup()
  res <- do.call(compose, .uf_args(s, detector = "builtin"))
  expect_identical(res$provenance$package, "Obligato")
  expect_true(all(c("classifier", "threshold", "feature_names", "type_median",
                    "phenotypes", "level_models") %in% names(res$detection_model)))
  expect_identical(res$detection_model$threshold, res$info$threshold)
})

test_that("invalid count and label inputs fail clearly", {
  s <- .uf_setup()
  bad <- s$query
  rownames(bad)[2] <- rownames(bad)[1]
  args_bad <- .uf_args(s, detector = "none")
  args_bad$query <- bad
  expect_error(do.call(compose, args_bad), "duplicated gene")

  args_labels <- .uf_args(s, detector = "none")
  args_labels$hc_labels <- LETTERS[1:3]
  expect_error(do.call(compose, args_labels), "hc_labels")
})
