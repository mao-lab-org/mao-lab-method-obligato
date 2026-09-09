################################################################################
# M1 GATE: the extraction must be behaviour-preserving.
#
# package_plan.md §10 Phase A: "lift functions verbatim ... Gate: existing outputs
# reproduce bit-identically. Any difference is a bug in the extraction, not an
# improvement." This test enforces that mechanically.
#
# It sources the ORIGINAL scripts into an isolated environment and compares each
# extracted function against its script counterpart on identical random inputs.
# When M1 is complete and the scripts call the package instead, this test becomes
# a regression guard against silent drift and should be kept, not deleted.
################################################################################

PROJ <- Sys.getenv("OBLIGATO_PROJECT_ROOT",
                   normalizePath(file.path("..", "..", ".."), mustWork = FALSE))

skip_if_no_scripts <- function() {
  f <- file.path(PROJ, "Phase0_Detection", "run_methodC_compose.R")
  skip_if_not(file.exists(f), "original scripts not reachable; set OBLIGATO_PROJECT_ROOT")
}

# Source only the function definitions from a script, without running its body.
# The scripts are top-to-bottom pipelines, so we parse and evaluate assignments only.
source_defs <- function(path) {
  e <- new.env(parent = globalenv())
  exprs <- parse(path)
  for (ex in exprs) {
    if (is.call(ex) && length(ex) >= 3 &&
        as.character(ex[[1]])[1] %in% c("<-", "=") &&
        is.call(ex[[3]]) && identical(as.character(ex[[3]][[1]])[1], "function")) {
      try(eval(ex, envir = e), silent = TRUE)
    }
  }
  e
}

test_that("logit-normal model functions are byte-identical to the scripts", {
  skip_if_no_scripts()
  skip_if_not_installed("corpcor")
  e <- source_defs(file.path(PROJ, "Phase0_Detection", "run_methodC_compose.R"))

  set.seed(1)
  X <- matrix(runif(40 * 6, -1, 1), nrow = 40)
  Y <- to_logit(X)

  expect_identical(to_logit(X),                       e$to_logit(X))
  expect_identical(pair_label("B naive", "CD4 TFH"),  e$pair_label("B naive", "CD4 TFH"))
  expect_identical(pair_label("CD4 TFH", "B naive"),  e$pair_label("CD4 TFH", "B naive"))

  a <- fit_ln_shrink(Y); b <- e$fit_ln_shrink(Y)
  expect_equal(a$mu,        b$mu,        tolerance = 0)
  expect_equal(a$Sigma_inv, b$Sigma_inv, tolerance = 0)
  expect_equal(a$log_det,   b$log_det,   tolerance = 0)

  expect_identical(quad_form_rows(Y, a$mu, a$Sigma_inv),
                   e$quad_form_rows(Y, b$mu, b$Sigma_inv))
  expect_identical(mvn_ll_rows(Y, a), e$mvn_ll_rows(Y, b))
})

test_that("candidate construction is identical to the script", {
  skip_if_no_scripts()
  e <- source_defs(file.path(PROJ, "Phase0_Detection", "run_methodC_compose.R"))
  top5 <- paste0("type", 1:5)
  expect_identical(make_candidates(top5), e$make_candidates(top5))
})

test_that("pair_label is order-invariant (property, not just parity)", {
  set.seed(3)
  a <- sample(letters, 50, TRUE); b <- sample(letters, 50, TRUE)
  expect_identical(pair_label(a, b), pair_label(b, a))
})

test_that("to_logit is monotone and finite on the score range", {
  x <- seq(-1, 1, length.out = 201)
  y <- to_logit(matrix(x, ncol = 1))
  expect_true(all(is.finite(y)))
  expect_true(all(diff(as.vector(y)) > 0))
})

test_that("high-confidence singlet helpers match the harness", {
  skip_if_no_scripts()
  h <- file.path(PROJ, "Phase0_Detection", "eval_harness.R")
  skip_if_not(file.exists(h))
  e <- source_defs(h)
  set.seed(4)
  p <- c(runif(400, 0, .2), runif(100, .6, .95))
  expect_identical(.clip01(c(-1, 0, .5, 1, 2)), e$.clip01(c(-1, 0, .5, 1, 2)))
  expect_identical(.logit(c(.1, .5, .9)),       e$.logit(c(.1, .5, .9)))
})
