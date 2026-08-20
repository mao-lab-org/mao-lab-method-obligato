################################################################################
# M1b GATE: ecdf_threshold() must reproduce the inline crossover computation
# from run_compose_kingetal.R:304-316 exactly.
#
# The inline block is not a function in the script, so it cannot be sourced and
# compared the way M1's verbatim extractions are. Instead the script's exact
# arithmetic is transcribed here as an independent oracle, and the extracted
# function is required to be identical to it. If the extraction ever drifts from
# the script's formula, this fails.
#
# package_plan.md §10 Phase B: this is "the function whose absence allowed the
# 0.5 error". End-to-end reproduction of the King threshold (0.711) is deferred
# to the M2b replication gate, which reruns the full pipeline.
################################################################################

# Independent transcription of run_compose_kingetal.R:304-316.
oracle_crossover <- function(real_scores, sim_scores, fallback = 0.5) {
  ecdf_real <- ecdf(real_scores)
  ecdf_sim  <- ecdf(sim_scores)
  crossover_fn <- function(t) ecdf_sim(t) - (1 - ecdf_real(t))
  lo <- min(c(real_scores, sim_scores))
  hi <- max(c(real_scores, sim_scores))
  if (crossover_fn(lo) * crossover_fn(hi) < 0) {
    uniroot(crossover_fn, c(lo, hi))$root
  } else {
    fallback
  }
}

test_that("ecdf_threshold is identical to the script's inline crossover", {
  set.seed(11)
  # real cells: mostly singlets (low scores) with a doublet tail;
  # simulated doublets: shifted high. This is the regime the crossover is for.
  real <- c(rbeta(4000, 2, 6), rbeta(400, 6, 2))
  sim  <- rbeta(2000, 6, 2)
  expect_identical(ecdf_threshold(real, sim), oracle_crossover(real, sim))
})

test_that("ecdf_threshold matches the oracle across many random regimes", {
  set.seed(12)
  for (i in 1:25) {
    n_real <- sample(200:3000, 1); n_sim <- sample(200:3000, 1)
    real <- rbeta(n_real, runif(1, 1, 8), runif(1, 1, 8))
    sim  <- rbeta(n_sim,  runif(1, 1, 8), runif(1, 1, 8))
    expect_identical(ecdf_threshold(real, sim), oracle_crossover(real, sim))
  }
})

test_that("the returned threshold satisfies the defining property (FNR_sim == flag_real)", {
  set.seed(13)
  real <- c(rbeta(4000, 2, 6), rbeta(400, 6, 2))
  sim  <- rbeta(2000, 6, 2)
  t_star <- ecdf_threshold(real, sim)
  fnr_sim   <- ecdf(sim)(t_star)          # sim doublets below t* (missed)
  flag_real <- 1 - ecdf(real)(t_star)     # real cells above t* (flagged)
  # Both ecdfs are step functions, so the crossover resolves the equality only to
  # within their jump size; bound the absolute gap by a couple of steps.
  expect_lt(abs(fnr_sim - flag_real), 2 / min(length(real), length(sim)))
  expect_gte(t_star, min(c(real, sim)))
  expect_lte(t_star, max(c(real, sim)))
})

test_that("degenerate input (no sign change) returns the fallback, and the fallback is settable", {
  const <- rep(0.5, 100)
  expect_identical(ecdf_threshold(const, const), 0.5)
  expect_identical(ecdf_threshold(const, const, fallback = 0.42), 0.42)
})
