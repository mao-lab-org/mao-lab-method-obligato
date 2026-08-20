##############################################################################
# Detection features and PhiSpace scoring.
# Extracted VERBATIM from run_methodC_compose.R (M1).
##############################################################################

# M1 DEVIATION — the only one, and it is not a behaviour change.
# In the scripts this function read `PHENO` from the enclosing script scope, where it
# was set from command-line arguments. A hidden global makes the function uncallable
# inside a package and is exactly the reproducibility defect packaging removes, so the
# dependency is promoted to an explicit argument. For any given `phenotypes` the
# computation is byte-identical to the script; only the source of that value changes.
# There is deliberately NO default: silently guessing the phenotype column is how
# `Cell.class` vs `predicted.celltype.l1` passed zero cells (package_plan.md §9).
score_phispace <- function(counts, ref_sce, phenotypes, label = "") {
  q <- SingleCellExperiment(assays = list(counts = counts))
  q <- scranTransf(q)
  q <- PhiSpace(reference = ref_sce, query = q, phenotypes = phenotypes,
                refAssay = "logcounts", queryAssay = "logcounts",
                regMethod = "PLS", center = TRUE, scale = FALSE,
                storeUnNorm = TRUE)
  S <- reducedDim(q, "PhiSpace")
  cat(sprintf("  [score%s] %d cells x %d types\n", label, nrow(S), ncol(S)))
  S
}

compute_features <- function(S, sing_models, pair_models) {
  Y <- to_logit(S)
  N <- nrow(S); Tn <- ncol(S)

  ord  <- t(apply(S, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
  top1 <- ord[, 1]; top2 <- ord[, 2]
  gap  <- top1 - top2
  Sx <- S - apply(S, 1, max)
  P  <- exp(Sx); P <- P / rowSums(P)
  ent <- -rowSums(P * log(P + 1e-12))
  eff_n <- exp(ent)

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

  data.frame(
    max_score   = top1,
    top_gap     = gap,
    entropy     = ent,
    eff_n       = eff_n,
    sing_LL     = max_sing_LL,
    dbl_LL      = max_dbl_LL,
    ll_diff     = max_dbl_LL - max_sing_LL,
    min_maha    = min_maha
  )
}

##############################################################################
# Data-driven detection threshold: the balanced ecdf crossover.
#
# Extracted from run_compose_kingetal.R:304-316 (identical block in
# run_compose_kingetal_full.R). package_plan.md §10 Phase B names this the
# highest-value extraction: as free-standing inline code the threshold could be
# bypassed by retyping a number (the `score > 0.5` error the package exists to
# prevent). As a function with the crossover as its only behaviour, it cannot.
#
# It finds t* where the false-negative rate on simulated doublets equals the
# flagging rate on real cells, i.e. ecdf_sim(t*) = 1 - ecdf_real(t*). This adapts
# to the actual score distributions without user input.
#
# Behaviour-preserving deviations from the inline block, both minor and documented:
#   * the hardcoded 0.5 fallback is exposed as the `fallback` argument, defaulting
#     to the same 0.5, so an override is recorded rather than edited into a copy;
#   * the script's cat() progress logging is dropped (it is a side effect, not part
#     of the returned value the parity test pins).
# For any given inputs the returned threshold is identical to the script's.
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

##############################################################################
# Train the xgboost detection classifier on singlet vs doublet features (M1b).
#
# Extracted from run_compose_kingetal.R:265-276 (Section 4), identical to
# run_methodC_compose.R's process_batch training block. Wraps the locked balanced
# 1:1 subsample and the locked xgboost hyperparameters.
#
# DETERMINISM — a sanctioned behaviour change, decided 2026-08-20 (package_plan
# §6, §10 Phase B). The scripts trained under xgboost's defaults: multithreaded
# and with no `seed`, which is NON-reproducible (two runs on identical inputs
# differ by up to ~0.27 in predicted score, moving the ecdf threshold and the
# flagged set). A package whose purpose is results that explain themselves cannot
# ship a non-deterministic default, so this function defaults to `nthread = 1` and
# a fixed `seed`, seeding BOTH the balancing subsample (R RNG) and xgboost's
# internal RNG. Consequence: it does not bit-reproduce the OLD stored detection
# scores (a non-deterministic draw); the M2b gate freezes a NEW deterministic
# reference and checks the old numbers only within tolerance. Composition, models
# and metrics are unaffected — they never depended on xgboost.
#
# Faithful-but-parameterised: with the defaults below the hyperparameters and the
# balancing rule are exactly the scripts'; only reproducibility is added.
train_detector <- function(F_sing, F_dbl,
                           params = list(objective = "binary:logistic",
                                         eval_metric = "auc",
                                         max_depth = 4, eta = 0.1,
                                         subsample = 0.8, colsample_bytree = 0.8),
                           nrounds = 200,
                           seed = 1L,
                           nthread = 1L) {
  # Balanced 1:1, bidirectional: subsample the larger class down to the smaller.
  set.seed(seed)
  n_bal <- min(nrow(F_sing), nrow(F_dbl))
  Fs <- F_sing[sample(nrow(F_sing), n_bal), , drop = FALSE]
  Fd <- F_dbl[sample(nrow(F_dbl),  n_bal), , drop = FALSE]

  X_train <- as.matrix(rbind(Fs, Fd))
  y_train <- c(rep(0L, nrow(Fs)), rep(1L, nrow(Fd)))

  p <- params
  p$nthread <- nthread
  p$seed    <- seed

  dtrain <- xgb.DMatrix(X_train, label = y_train)
  xgb.train(params = p, data = dtrain, nrounds = nrounds, verbose = 0)
}
