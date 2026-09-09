##############################################################################
# Detection metrics.
# Extracted VERBATIM from eval_harness.R (M1).
##############################################################################

# Trapezoidal AUC.
.trap_auc <- function(x, y) {
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}

# Compute precision-recall curve and AUPRC.
compute_auprc <- function(truth, scores) {
  if (length(unique(truth)) < 2) return(list(auprc = NA_real_, curve = NULL))
  ord <- order(scores, decreasing = TRUE)
  t <- truth[ord]
  tp <- cumsum(t == 1)
  fp <- cumsum(t == 0)
  P  <- sum(truth == 1)
  precision <- tp / (tp + fp)
  recall    <- tp / P
  # AUPRC via trapezoidal on (recall, precision) sorted by recall
  auprc <- .trap_auc(c(0, recall), c(precision[1], precision))
  list(auprc = auprc,
       curve = data.frame(recall = recall, precision = precision))
}

##############################################################################
# Composition metrics (M1b).
#
# Phase B "wrap, don't copy": these are the paper's claimed evaluation
# contributions (paper_plan §1.4e) with no single reusable definition in the
# scripts, so they are written to clean signatures and pinned to the exact
# hand-computed worked example in Phase0_Detection/scdblfinder_composition_plan.md
# §8 (5 doublets; note package_plan §8.1's "6" is a typo — the numbers 0.60 /
# 0.40 / 0.667 are the 5-doublet values). The decomposition identity in §8.3 is
# asserted by the test, per the instruction to port it as the metrics test.
##############################################################################

# Exact-pair composition accuracy ("Correct-L1"): a call is correct when the
# unordered predicted type pair equals the unordered true type pair. Returns both
# the realistic (overall) and balanced (equal weight per true pair type) values —
# the locked reporting rule — plus per-stratum accuracy when `lineage` is given.
composition_accuracy <- function(pred_A, pred_B, true_A, true_B, lineage = NULL) {
  n <- length(true_A)
  stopifnot(length(pred_A) == n, length(pred_B) == n, length(true_B) == n)
  correct <- mapply(function(pa, pb, ta, tb) setequal(c(pa, pb), c(ta, tb)),
                    pred_A, pred_B, true_A, true_B, USE.NAMES = FALSE)

  pair_key <- pair_label(true_A, true_B)          # unordered true-pair key
  balanced <- mean(tapply(correct, pair_key, mean))

  out <- list(
    correct  = correct,
    overall  = mean(correct),
    balanced = balanced,
    n        = n)
  if (!is.null(lineage)) out$by_stratum <- tapply(correct, lineage, mean)
  out
}

# Majority decode of a grouping: the majority SINGLET cell type in each group.
# Reproduces the eval harness's `phi` exactly (table + max.col, first-tie). This
# is the reproducible bridge from a cluster-output method's groups to L1 types.
majority_decode <- function(grouping, types) {
  ok  <- !is.na(grouping) & !is.na(types)
  tab <- table(grouping[ok], types[ok])
  stats::setNames(colnames(tab)[max.col(tab, ties.method = "first")], rownames(tab))
}

# Maximum achievable accuracy: the best exact-pair accuracy any composition rule
# could reach on a given grouping, independent of predictions. Two definitions
# (both reported, per scdblfinder_composition_plan §10.6):
#   "vocabulary"  (PRIMARY, the true upper bound): both true types are the decode
#                 of SOME group, so the answer space contains the right answer.
#   "constituent" (tighter): the two constituents lie in different groups AND
#                 decoding those groups reproduces the truth. Needs group_A/group_B
#                 and is reported over the both-constituents-present subset.
# `decode` is a named vector group -> type (e.g. from majority_decode()).
max_achievable_accuracy <- function(true_A, true_B, decode,
                                    group_A = NULL, group_B = NULL,
                                    lineage = NULL,
                                    definition = c("vocabulary", "constituent")) {
  definition <- match.arg(definition)

  if (definition == "vocabulary") {
    reach <- unique(unname(decode))
    ma    <- (true_A %in% reach) & (true_B %in% reach)
    keep  <- rep(TRUE, length(ma))
  } else {
    if (is.null(group_A) || is.null(group_B)) {
      stop("definition = \"constituent\" requires group_A and group_B.", call. = FALSE)
    }
    both       <- !is.na(group_A) & !is.na(group_B)
    resolvable <- both & (group_A != group_B)
    dA <- decode[group_A]; dB <- decode[group_B]
    decodes_ok <- mapply(function(a, b, x, y)
                           !is.na(a) && !is.na(b) && setequal(c(a, b), c(x, y)),
                         dA, dB, true_A, true_B, USE.NAMES = FALSE)
    ma   <- resolvable & decodes_ok
    keep <- both
  }

  out <- list(max_achievable = mean(ma[keep]), definition = definition, n = sum(keep))
  if (!is.null(lineage)) out$by_stratum <- tapply(ma[keep], lineage[keep], mean)
  out
}
