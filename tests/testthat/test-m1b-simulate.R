################################################################################
# M1b GATE: simulate_training_doublets() reproduces Section 1 of
# run_compose_kingetal.R:184-201, given the same seed (determinism is the
# sanctioned change; see the function header).
################################################################################

oracle_sim <- function(counts, hc_idx, hc_types, types, n_per_pair, seed) {
  set.seed(seed)
  type_to_cells <- split(seq_along(hc_idx), hc_types)
  het_pairs <- combn(types, 2, simplify = FALSE)
  dblA <- integer(0); dblB <- integer(0); dbl_pl <- character(0)
  for (pr in het_pairs) {
    a <- pr[1]; b <- pr[2]
    if (length(type_to_cells[[a]]) < 1 || length(type_to_cells[[b]]) < 1) next
    ca <- hc_idx[sample(type_to_cells[[a]], n_per_pair, replace = TRUE)]
    cb <- hc_idx[sample(type_to_cells[[b]], n_per_pair, replace = TRUE)]
    dblA <- c(dblA, ca); dblB <- c(dblB, cb)
    dbl_pl <- c(dbl_pl, rep(pair_label(a, b), n_per_pair))
  }
  list(A = dblA, B = dblB, pl = dbl_pl)
}

make_counts <- function(n_genes, n_cells, seed = 5) {
  set.seed(seed)
  matrix(rpois(n_genes * n_cells, 3), nrow = n_genes,
         dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells))))
}

test_that("simulate_training_doublets matches the script's Section 1 (same seed)", {
  types <- c("A", "B", "C", "D")
  counts <- make_counts(30, 240)
  hc_idx <- 1:240
  hc_types <- rep(types, each = 60)
  sim <- simulate_training_doublets(counts, hc_idx, hc_types, types, n_per_pair = 40, seed = 3L)
  orc <- oracle_sim(counts, hc_idx, hc_types, types, 40, 3L)

  expect_identical(sim$idx_A, orc$A)
  expect_identical(sim$idx_B, orc$B)
  expect_identical(sim$pair_label, orc$pl)
  # counts are the exact element-wise sum of the two constituents
  expect_equal(unname(sim$counts), unname(counts[, orc$A] + counts[, orc$B]))
  # C(4,2)=6 pairs x 40 doublets
  expect_equal(ncol(sim$counts), 6 * 40)
  expect_identical(sort(unique(sim$pair_label)),
                   sort(vapply(combn(types, 2, simplify = FALSE),
                               function(p) pair_label(p[1], p[2]), character(1))))
})

test_that("simulate_training_doublets is deterministic and seed-sensitive", {
  types <- c("A", "B", "C")
  counts <- make_counts(20, 150)
  hc_idx <- 1:150; hc_types <- rep(types, each = 50)
  s1 <- simulate_training_doublets(counts, hc_idx, hc_types, types, 30, seed = 1L)
  s2 <- simulate_training_doublets(counts, hc_idx, hc_types, types, 30, seed = 1L)
  s3 <- simulate_training_doublets(counts, hc_idx, hc_types, types, 30, seed = 2L)
  expect_identical(s1$idx_A, s2$idx_A)
  expect_false(identical(s1$idx_A, s3$idx_A))
})

test_that("a type with no high-confidence singlets is skipped", {
  types <- c("A", "B", "C")     # C present in `types` but absent from hc_types
  counts <- make_counts(20, 100)
  hc_idx <- 1:100; hc_types <- rep(c("A", "B"), each = 50)
  sim <- simulate_training_doublets(counts, hc_idx, hc_types, types, 30, seed = 1L)
  expect_identical(unique(sim$pair_label), pair_label("A", "B"))
})
