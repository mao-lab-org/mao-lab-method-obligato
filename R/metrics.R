# Detection and composition evaluation metrics.

# Trapezoidal AUC.
.trap_auc <- function(x, y) {
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}

#' Compute area under a precision-recall curve
#'
#' Tied scores are evaluated as one threshold, making the result independent of
#' the input row order. The area is calculated by trapezoidal integration.
#' @param truth Logical or binary vector, where one denotes a doublet.
#' @param scores Numeric detection scores; larger values indicate doublets.
#' @return A list with `auprc` and a threshold-level `curve`.
#' @export
compute_auprc <- function(truth, scores) {
  if (length(truth) != length(scores) || !length(truth) ||
      anyNA(truth) || anyNA(scores) || any(!is.finite(scores)) ||
      any(!truth %in% c(FALSE, TRUE, 0, 1))) {
    stop("`truth` and `scores` must be equal-length, complete binary/numeric vectors.",
         call. = FALSE)
  }
  truth <- as.integer(truth)
  if (length(unique(truth)) < 2L) {
    return(list(auprc = NA_real_, curve = NULL))
  }
  ord <- order(scores, decreasing = TRUE)
  s <- scores[ord]
  t <- truth[ord]
  tp <- cumsum(t == 1L)
  fp <- cumsum(t == 0L)
  keep <- !duplicated(s, fromLast = TRUE)
  tp <- tp[keep]; fp <- fp[keep]; thresholds <- s[keep]
  precision <- tp / (tp + fp)
  recall <- tp / sum(truth == 1L)
  auprc <- .trap_auc(c(0, recall), c(precision[1L], precision))
  list(auprc = auprc,
       curve = data.frame(threshold = thresholds, recall = recall,
                          precision = precision))
}

# Exact-pair metrics treat constituent order as irrelevant.

#' Evaluate exact unordered-pair composition accuracy
#'
#' Missing predictions are counted as incorrect. Balanced accuracy gives each
#' observed true pair equal weight, irrespective of its number of cells.
#' @param pred_A,pred_B Predicted constituent labels.
#' @param true_A,true_B Ground-truth constituent labels.
#' @param lineage Optional stratum label with one value per doublet.
#' @return A list containing per-cell correctness, overall and pair-balanced
#'   accuracy, sample size, and optional stratum accuracies.
#' @export
composition_accuracy <- function(pred_A, pred_B, true_A, true_B, lineage = NULL) {
  n <- length(true_A)
  pred_A <- as.character(pred_A); pred_B <- as.character(pred_B)
  true_A <- as.character(true_A); true_B <- as.character(true_B)
  if (!n || any(lengths(list(pred_A, pred_B, true_B)) != n) ||
      anyNA(true_A) || anyNA(true_B) ||
      (!is.null(lineage) && length(lineage) != n)) {
    stop("Composition vectors must be non-empty and have matching lengths; truth cannot be missing.",
         call. = FALSE)
  }
  correct <- mapply(function(pa, pb, ta, tb) {
    if (is.na(pa) || is.na(pb)) FALSE else setequal(c(pa, pb), c(ta, tb))
  }, pred_A, pred_B, true_A, true_B, USE.NAMES = FALSE)

  pair_key <- pair_label(true_A, true_B)
  balanced <- mean(tapply(correct, pair_key, mean))
  out <- list(correct = correct, overall = mean(correct),
              balanced = balanced, n = n)
  if (!is.null(lineage)) out$by_stratum <- tapply(correct, lineage, mean)
  out
}

#' Decode groups by their majority singlet type
#' @param grouping Group identifier for each singlet.
#' @param types Cell-type label for each singlet.
#' @return A named character vector mapping groups to cell types.
#' @export
majority_decode <- function(grouping, types) {
  if (length(grouping) != length(types)) {
    stop("`grouping` and `types` must have equal length.", call. = FALSE)
  }
  ok <- !is.na(grouping) & !is.na(types)
  if (!any(ok)) return(stats::setNames(character(), character()))
  tab <- table(grouping[ok], types[ok])
  stats::setNames(colnames(tab)[max.col(tab, ties.method = "first")], rownames(tab))
}

#' Calculate the maximum accuracy permitted by a grouping
#'
#' The `"vocabulary"` definition asks whether both true types occur in the
#' decoded group vocabulary. The `"constituent"` definition additionally
#' requires the two observed constituents to occupy different groups that decode
#' to the true unordered pair.
#' @param true_A,true_B Ground-truth constituent labels.
#' @param decode Named character vector mapping group identifiers to cell types.
#' @param group_A,group_B Group identifiers for the two constituents; required for
#'   `definition = "constituent"`.
#' @param lineage Optional stratum label.
#' @param definition Ceiling definition to calculate.
#' @return A list containing the ceiling, definition, denominator, and optional
#'   stratum-specific ceilings.
#' @export
max_achievable_accuracy <- function(true_A, true_B, decode,
                                    group_A = NULL, group_B = NULL,
                                    lineage = NULL,
                                    definition = c("vocabulary", "constituent")) {
  definition <- match.arg(definition)
  n <- length(true_A)
  if (!n || length(true_B) != n || anyNA(true_A) || anyNA(true_B) ||
      is.null(names(decode)) || anyNA(names(decode)) || anyDuplicated(names(decode)) ||
      (!is.null(lineage) && length(lineage) != n)) {
    stop("Truth vectors must be complete and matched; `decode` must be uniquely named.",
         call. = FALSE)
  }
  if (definition == "constituent" && (is.null(group_A) || is.null(group_B))) {
    stop("definition = \"constituent\" requires group_A and group_B.",
         call. = FALSE)
  }
  if (definition == "constituent" &&
      (length(group_A) != n || length(group_B) != n)) {
    stop("`group_A` and `group_B` must have one value per true pair.",
         call. = FALSE)
  }
  true_A <- as.character(true_A); true_B <- as.character(true_B)
  decode <- stats::setNames(as.character(decode), names(decode))
  if (!is.null(group_A)) group_A <- as.character(group_A)
  if (!is.null(group_B)) group_B <- as.character(group_B)

  if (definition == "vocabulary") {
    reach <- unique(unname(decode))
    ma    <- (true_A %in% reach) & (true_B %in% reach)
    keep  <- rep(TRUE, length(ma))
  } else {
    both       <- !is.na(group_A) & !is.na(group_B)
    resolvable <- both & (group_A != group_B)
    dA <- decode[group_A]; dB <- decode[group_B]
    decodes_ok <- mapply(function(a, b, x, y)
                           !is.na(a) && !is.na(b) && setequal(c(a, b), c(x, y)),
                         dA, dB, true_A, true_B, USE.NAMES = FALSE)
    ma   <- resolvable & decodes_ok
    keep <- both
  }

  value <- if (any(keep)) mean(ma[keep]) else NA_real_
  out <- list(max_achievable = value, definition = definition, n = sum(keep))
  if (!is.null(lineage)) {
    out$by_stratum <- if (any(keep)) tapply(ma[keep], lineage[keep], mean) else NULL
  }
  out
}
