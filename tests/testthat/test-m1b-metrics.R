################################################################################
# M1b GATE: composition metrics against the hand-computed worked example in
# Phase0_Detection/scdblfinder_composition_plan.md §8. Every number is exact and
# verifiable by hand; the functions must reproduce them, and the loss
# decomposition (§8.3) must balance.
################################################################################

# --- The worked example -------------------------------------------------------
# Four types; a cluster method with c1 impure (60% CD4T / 40% CD8T), c2=Mono,
# c3=B, c4=B (c3/c4 redundant).
sing_group <- c(rep("c1", 10), rep("c2", 5), rep("c3", 5), rep("c4", 5))
sing_type  <- c(rep("CD4T", 6), rep("CD8T", 4), rep("Mono", 5), rep("B", 5), rep("B", 5))

# Five ground-truth heterotypic doublets (§8.1).
true_A <- c("CD4T", "CD8T", "CD4T", "B",    "B")
true_B <- c("Mono", "Mono", "CD8T", "Mono", "CD4T")
grp_A  <- c("c1",   "c1",   "c1",   "c3",   "c3")
grp_B  <- c("c2",   "c2",   "c1",   "c2",   "c1")
lineage <- c("cross", "cross", "within", "cross", "cross")

# Decode of the method's predictions to L1 (§8.2 "decode of prediction" column).
pred_A <- c("CD4T", "CD4T", "CD4T", "B",    "Mono")
pred_B <- c("Mono", "Mono", "Mono", "Mono", "B")

test_that("majority_decode reproduces the harness phi (first-tie majority)", {
  dec <- majority_decode(sing_group, sing_type)
  expect_identical(unname(dec[c("c1", "c2", "c3", "c4")]), c("CD4T", "Mono", "B", "B"))
})

test_that("composition_accuracy reproduces the worked example (§8.3, §8.4)", {
  ca <- composition_accuracy(pred_A, pred_B, true_A, true_B, lineage = lineage)
  expect_equal(ca$overall, 0.40)                       # 2/5, headline Correct-L1
  expect_equal(ca$balanced, 0.40)                      # every true pair unique
  expect_equal(unname(ca$by_stratum["within"]), 0.00)  # {d3}
  expect_equal(unname(ca$by_stratum["cross"]), 0.50)   # {d1,d2,d4,d5} -> 2/4
  expect_identical(ca$correct, c(TRUE, FALSE, FALSE, TRUE, FALSE))
})

test_that("max_achievable_accuracy = 0.60 under BOTH definitions (§8.3)", {
  dec <- majority_decode(sing_group, sing_type)
  mv <- max_achievable_accuracy(true_A, true_B, dec, lineage = lineage,
                                definition = "vocabulary")
  mc <- max_achievable_accuracy(true_A, true_B, dec, grp_A, grp_B, lineage = lineage,
                                definition = "constituent")
  expect_equal(mv$max_achievable, 0.60)
  expect_equal(mc$max_achievable, 0.60)
  # vocabulary strata: within {d3}=0, cross {d1,d2,d4,d5}={T,F,T,T}=3/4
  expect_equal(unname(mv$by_stratum["within"]), 0.00)
  expect_equal(unname(mv$by_stratum["cross"]),  0.75)
  expect_equal(unname(mc$by_stratum["within"]), 0.00)
  expect_equal(unname(mc$by_stratum["cross"]),  0.75)
})

test_that("the loss decomposition balances (§8.3)", {
  dec <- majority_decode(sing_group, sing_type)
  max_a  <- max_achievable_accuracy(true_A, true_B, dec, definition = "vocabulary")$max_achievable
  corr_l1 <- composition_accuracy(pred_A, pred_B, true_A, true_B)$overall
  representational_loss <- 1 - max_a          # 0.40
  identification_loss   <- max_a - corr_l1    # 0.20
  expect_equal(representational_loss, 0.40)
  expect_equal(identification_loss, 0.20)
  expect_equal(representational_loss + identification_loss, 1 - corr_l1)  # balances
})

test_that("constituent definition requires groups; vocabulary does not", {
  dec <- majority_decode(sing_group, sing_type)
  expect_error(max_achievable_accuracy(true_A, true_B, dec, definition = "constituent"),
               "requires group_A and group_B")
  expect_silent(max_achievable_accuracy(true_A, true_B, dec, definition = "vocabulary"))
})

test_that("reference-based arm: L1 grouping makes every heterotypic pair achievable", {
  # DblLik: groups ARE L1 labels, decode is identity -> reach is all types (§8.5).
  dec_identity <- setNames(c("CD4T", "CD8T", "Mono", "B"), c("CD4T", "CD8T", "Mono", "B"))
  mv <- max_achievable_accuracy(true_A, true_B, dec_identity, definition = "vocabulary")
  expect_equal(mv$max_achievable, 1.00)
})
