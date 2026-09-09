#' Print an Obligato composition result
#' @param x An `obligato_composition` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.obligato_composition <- function(x, ...) {
  cat("<obligato_composition>\n")
  cat("  cells:       ", x$info$n_cells, "\n", sep = "")
  cat("  composition: ", x$info$output, " (", x$info$candidates, ")\n", sep = "")
  cat("  detector:    ", x$info$detector, "\n", sep = "")
  if (!is.null(x$flag)) {
    cat("  flagged:     ", sum(x$flag), " (",
        sprintf("%.1f%%", 100 * mean(x$flag)), ")\n", sep = "")
  }
  invisible(x)
}

#' Print an Obligato reference-free result
#' @param x An `obligato_reffree` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.obligato_reffree <- function(x, ...) {
  cat("<obligato_reffree>\n")
  cat("  cells:       ", x$info$n_cells, "\n", sep = "")
  cat("  rounds:      ", x$info$n_rounds + 1L, "\n", sep = "")
  cat("  clusters:    ", x$info$n_final_clusters, "\n", sep = "")
  cat("  flagged:     ", sum(x$flag), " (",
      sprintf("%.1f%%", 100 * mean(x$flag)), ")\n", sep = "")
  invisible(x)
}
