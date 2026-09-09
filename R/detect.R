# Internal PhiSpace scoring and detector features.
#
# The phenotype column is always explicit: guessing it can silently select the
# wrong reference annotation.
score_phispace <- function(counts, ref_sce, phenotypes, label = "") {
  q <- SingleCellExperiment(assays = list(counts = counts))
  q <- scranTransf(q)
  q <- PhiSpace(reference = ref_sce, query = q, phenotypes = phenotypes,
                refAssay = "logcounts", queryAssay = "logcounts",
                regMethod = "PLS", center = TRUE, scale = FALSE,
                storeUnNorm = TRUE)
  S <- reducedDim(q, "PhiSpace")
  S
}

# Library-size features restore count-depth information removed by PhiSpace
# normalization: log total counts, detected genes, and log total counts relative
# to the median of the cell's highest-scoring type. They affect detection only.
libsize_features <- function(counts, top1_type, type_median = NULL) {
  lt <- log1p(Matrix::colSums(counts))
  ng <- as.numeric(Matrix::colSums(counts > 0))
  out <- data.frame(log_total = lt, n_genes = ng)
  if (!is.null(type_median)) {
    tr <- lt - type_median[top1_type]
    tr[is.na(tr)] <- 0
    out$typerel <- tr
  }
  out
}

# Assign each cell to its highest-scoring type. Keeping this in one helper makes
# the type-relative library-size definition identical in training and inference.
top1_type <- function(S) {
  if (is.null(colnames(S)) || ncol(S) < 1L) {
    stop("`S` must have at least one named score column.", call. = FALSE)
  }
  colnames(S)[max.col(S, ties.method = "first")]
}

# Fit the type-relative library-size baseline on the score-derived top-1 types,
# matching the validated A2 ablation. Do not group by supplied hc labels: those
# labels can differ from the jointly normalised training scores used by detection.
libsize_type_medians <- function(counts, S) {
  if (ncol(counts) != nrow(S)) {
    stop("`counts` columns must match `S` rows.", call. = FALSE)
  }
  base::tapply(log1p(Matrix::colSums(counts)), top1_type(S), stats::median)
}

compute_features <- function(S, sing_models, pair_models, suffix = "") {
  Y <- to_logit(S)
  N <- nrow(S); Tn <- ncol(S)

  ord  <- t(apply(S, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
  top1 <- ord[, 1]; top2 <- ord[, 2]
  gap  <- top1 - top2
  Sx <- S - apply(S, 1, max)
  P  <- exp(Sx); P <- P / rowSums(P)
  ent <- -rowSums(P * log(P + 1e-12))

  sing_LL  <- matrix(-Inf, N, length(sing_models))
  sing_q   <- matrix( Inf, N, length(sing_models))
  for (j in seq_along(sing_models)) {
    m <- sing_models[[j]]
    sing_q[, j]  <- quad_form_rows(Y, m$mu, m$Sigma_inv)
    sing_LL[, j] <- -0.5 * (Tn * log(2 * pi) + m$log_det + sing_q[, j])
  }
  max_sing_LL <- apply(sing_LL, 1, max)
  min_maha    <- sqrt(apply(sing_q, 1, min))

  dbl_LL <- matrix(-Inf, N, length(pair_models))
  for (j in seq_along(pair_models)) dbl_LL[, j] <- mvn_ll_rows(Y, pair_models[[j]])
  max_dbl_LL <- apply(dbl_LL, 1, max)

  feats <- data.frame(
    max_score   = top1,
    top_gap     = gap,
    entropy     = ent,
    sing_LL     = max_sing_LL,
    dbl_LL      = max_dbl_LL,
    ll_diff     = max_dbl_LL - max_sing_LL,
    min_maha    = min_maha
  )
  # Per-level suffix (e.g. "_l1", "_l2") so features from concatenated annotation
  # levels stay distinct. Library size is added once at the compose() level (it is
  # level-independent), not here.
  if (nzchar(suffix)) names(feats) <- paste0(names(feats), suffix)
  feats
}

# Balanced ECDF crossover threshold. It finds the score where the false-negative
# rate among simulated doublets equals the flagging rate among query cells.
ecdf_threshold <- function(real_scores, sim_scores, fallback = 0.5) {
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

# Train the xgboost detector on a balanced 1:1 singlet/doublet sample. A
# single thread and explicit seeds provide reproducible package behavior.
train_detector <- function(F_sing, F_dbl,
                           params = list(objective = "binary:logistic",
                                         eval_metric = "auc",
                                         max_depth = 4, eta = 0.1,
                                         subsample = 0.8, colsample_bytree = 0.8),
                           nrounds = 200,
                           seed = 1L,
                           nthread = 1L) {
  if (!is.data.frame(F_sing) || !is.data.frame(F_dbl) ||
      nrow(F_sing) < 1L || nrow(F_dbl) < 1L ||
      !identical(names(F_sing), names(F_dbl))) {
    stop("Training feature tables must be non-empty and have identical columns.",
         call. = FALSE)
  }
  nrounds <- .validate_scalar_integer(nrounds, "nrounds", 1L)
  nthread <- .validate_scalar_integer(nthread, "nthread", 1L)
  seed <- .validate_scalar_integer(seed, "seed", 0L)

  n_bal <- min(nrow(F_sing), nrow(F_dbl))
  sampled <- .with_seed(seed, list(
    sing = sample(nrow(F_sing), n_bal),
    dbl = sample(nrow(F_dbl), n_bal)))
  Fs <- F_sing[sampled$sing, , drop = FALSE]
  Fd <- F_dbl[sampled$dbl, , drop = FALSE]

  X_train <- as.matrix(rbind(Fs, Fd))
  y_train <- c(rep(0L, nrow(Fs)), rep(1L, nrow(Fd)))
  p <- params
  p$nthread <- nthread
  p$seed <- seed

  dtrain <- xgb.DMatrix(X_train, label = y_train)
  xgb.train(params = p, data = dtrain, nrounds = nrounds, verbose = 0)
}
