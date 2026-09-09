# Internal validation and reproducibility helpers.

.validate_counts <- function(x, name = "query") {
  if (length(dim(x)) != 2L || nrow(x) < 1L || ncol(x) < 1L) {
    stop(sprintf("`%s` must be a non-empty genes-by-cells matrix-like object.", name),
         call. = FALSE)
  }
  vals <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
  if (!is.numeric(vals) || anyNA(vals) || any(!is.finite(vals)) || any(vals < 0)) {
    stop(sprintf("`%s` must contain finite, non-negative numeric counts.", name),
         call. = FALSE)
  }
  if (!is.null(rownames(x)) && anyDuplicated(rownames(x))) {
    stop(sprintf("`%s` has duplicated gene names; make row names unique first.", name),
         call. = FALSE)
  }
  invisible(x)
}

.validate_scores <- function(S, n, label = "score matrix", expected = NULL) {
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) != n || ncol(S) < 2L) {
    stop(sprintf("%s must be a numeric %d-row matrix with at least two cell-type columns.",
                 label, n), call. = FALSE)
  }
  if (is.null(colnames(S)) || anyNA(colnames(S)) || anyDuplicated(colnames(S))) {
    stop(sprintf("%s must have unique, non-missing cell-type column names.", label),
         call. = FALSE)
  }
  if (anyNA(S) || any(!is.finite(S))) {
    stop(sprintf("%s contains missing or non-finite values.", label), call. = FALSE)
  }
  if (!is.null(expected) && !identical(colnames(S), expected)) {
    stop(sprintf("%s has cell-type columns inconsistent with the primary score matrix.", label),
         call. = FALSE)
  }
  S
}

.validate_scalar_integer <- function(x, name, lower = 0L) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x != as.integer(x) || x < lower) {
    stop(sprintf("`%s` must be a single integer >= %d.", name, lower), call. = FALSE)
  }
  as.integer(x)
}

# Evaluate an expression under a local deterministic R seed without changing the
# caller's global RNG state. This also makes injectable scoring functions obey the
# public `seed` contract when their internals use randomised normalisation.
.with_seed <- function(seed, code) {
  seed <- .validate_scalar_integer(seed, "seed", 0L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.package_provenance <- function() {
  version <- tryCatch(as.character(utils::packageVersion("Obligato")),
                      error = function(e) NA_character_)
  list(package = "Obligato", version = version,
       generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
}
