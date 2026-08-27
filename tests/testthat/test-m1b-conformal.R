################################################################################
# M1b: split-conformal prediction sets for composition (S4 / Item 3).
# Unit tests for the pure functions (conformal_calibrate / conformal_sets /
# conformal_pairs) and an end-to-end compose(output = "conformal") smoke test
# with an injected deterministic scorer. The reproduction of the reported set
# sizes (B2 6.7 / B9 7.6 / PBMC 2.7) is a data-dependent integration gate, not a
# unit test.
################################################################################

# Build a P-pair log-likelihood matrix with a controllable true pair per row.
make_logh <- function(n, P, seed, true_mean = 3) {
  set.seed(seed)
  pm  <- paste0("p", seq_len(P))
  M   <- matrix(rnorm(n * P), n, P, dimnames = list(NULL, pm))
  tp  <- sample(pm, n, replace = TRUE)            # true pair per row
  for (i in seq_len(n)) M[i, tp[i]] <- rnorm(1, true_mean)  # true pair scores higher
  list(log_h = M, true_pairs = tp, pm = pm)
}

test_that("conformal_calibrate: q_hat increases as coverage increases (alpha down)", {
  d <- make_logh(500, 8, seed = 1)
  q <- conformal_calibrate(d$log_h, d$true_pairs, alpha = c(0.20, 0.10, 0.01))
  expect_length(q, 3)
  expect_true(q[["alpha_0.2"]] <= q[["alpha_0.1"]])
  expect_true(q[["alpha_0.1"]] <= q[["alpha_0.01"]])
})

test_that("conformal_sets: membership is loglik >= -q_hat and sizes match", {
  d <- make_logh(50, 6, seed = 2)
  q <- conformal_calibrate(d$log_h, d$true_pairs, alpha = 0.10)[[1]]
  cs <- conformal_sets(d$log_h, q)
  expect_length(cs$sets, 50)
  expect_equal(cs$sizes, lengths(cs$sets))
  # every member has loglik >= -q_hat; every non-member is below
  for (i in c(1, 17, 50)) {
    in_p  <- cs$sets[[i]]
    expect_true(all(d$log_h[i, in_p] >= -q))
    out_p <- setdiff(d$pm, in_p)
    if (length(out_p)) expect_true(all(d$log_h[i, out_p] < -q))
  }
})

test_that("split-conformal attains marginal coverage >= 1 - alpha (finite-sample band)", {
  # Calibration and query drawn exchangeably -> coverage guarantee holds.
  alpha <- 0.10
  d_cal <- make_logh(2000, 10, seed = 10)
  d_qry <- make_logh(2000, 10, seed = 11)
  cf <- conformal_pairs(d_cal$log_h, d_cal$true_pairs, d_qry$log_h, alpha = alpha)
  covered <- vapply(seq_along(cf$sets), function(i)
    d_qry$true_pairs[i] %in% cf$sets[[i]], logical(1))
  emp <- mean(covered)
  # marginal coverage in [1-alpha - slack, 1-alpha + slack]; conformal both-sides
  expect_gt(emp, 1 - alpha - 0.03)
  expect_lt(emp, 1 - alpha + 0.05)
})

test_that("true pair absent from models -> Inf nonconformity, conservative sets", {
  d <- make_logh(100, 5, seed = 3)
  tp_missing <- rep("not_a_pair", 100)             # no true pair is ever modelled
  q <- conformal_calibrate(d$log_h, tp_missing, alpha = 0.10)[[1]]
  expect_true(is.infinite(q))                       # q_hat = Inf
  cs <- conformal_sets(d$log_h, q)                  # every pair included
  expect_true(all(cs$sizes == 5))
})

test_that("conformal_calibrate validates inputs", {
  d <- make_logh(10, 4, seed = 4)
  expect_error(conformal_calibrate(d$log_h, d$true_pairs, alpha = 0),   "alpha")
  expect_error(conformal_calibrate(d$log_h, d$true_pairs, alpha = 1.2), "alpha")
  expect_error(conformal_calibrate(d$log_h, d$true_pairs[1:5]))          # length mismatch
})

# ---- end-to-end compose(output = "conformal") --------------------------------
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
  hc_labels <- rep(types, each = 20)
  list(query = query, reference = reference, hc = hc, hc_labels = hc_labels)
}

test_that("compose(output = 'conformal') returns sets, sizes, and q_hat", {
  s <- setup_query()
  res <- compose(query = s$query, reference = s$reference, phenotypes = "celltype",
                 hc_singlets = s$hc, hc_labels = s$hc_labels, n_per_pair = 30,
                 detector = "none", output = "conformal", alpha = 0.10,
                 score_fn = fake_scorer)
  expect_s3_class(res, "obligato_composition")
  expect_false(is.null(res$conformal))
  expect_equal(res$info$output, "conformal")
  expect_equal(res$info$alpha, 0.10)
  # one size per query cell; composition is long-format members of each set
  expect_length(res$conformal$sizes, 200)
  expect_true(all(c("cell", "c1", "c2", "pair", "ll") %in% names(res$composition)))
  # rows == total membership across cells
  expect_equal(nrow(res$composition), sum(res$conformal$sizes))
  # every listed pair's ll clears the -q_hat bar
  expect_true(all(res$composition$ll >= -res$conformal$q_hat))
})

test_that("compose conformal is deterministic given seed", {
  s <- setup_query()
  args <- list(query = s$query, reference = s$reference, phenotypes = "celltype",
               hc_singlets = s$hc, hc_labels = s$hc_labels, n_per_pair = 30,
               detector = "none", output = "conformal", alpha = 0.10,
               score_fn = fake_scorer, seed = 7L)
  r1 <- do.call(compose, args); r2 <- do.call(compose, args)
  expect_identical(r1$conformal$q_hat, r2$conformal$q_hat)
  expect_identical(r1$conformal$sizes, r2$conformal$sizes)
})
