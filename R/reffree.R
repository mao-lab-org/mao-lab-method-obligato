##############################################################################
# Reference-free mode (M2c): the iterative detect -> remove -> recluster loop.
#
# Wraps the procedure in Phase0_Detection/run_methodC_reffree_compose.R so the
# reported reference-free results are reproducible from the package. No external
# labelled reference: each round clusters the current doublet-free pool, treats the
# cluster IDs as "cell types", uses the cluster-labelled pool as the PhiSpace
# reference, simulates doublets between clusters, fits the same logit-normal models,
# trains the same detector, flags at the DATA-DRIVEN ecdf-crossover threshold (NOT a
# fixed 0.5; decided 2026-08-24), removes, and repeats. Composition afterwards ranks
# every fitted cluster-pair model by likelihood (all_pair_ll), as in the script.
#
# It reuses the packaged building blocks (to_logit, fit_ln_shrink, compute_features,
# train_detector, ecdf_threshold, all_pair_ll); only the clustering and the
# cluster-pair simulation are reference-free-specific. Clustering is INJECTABLE
# (`cluster_fn`) so the package does not hard-depend on Seurat — the default uses
# Seurat when available; tests inject a deterministic stub.
##############################################################################

# Default clustering: Seurat PCA -> SNN -> Louvain. Returns a character cluster id
# per cell (columns of `counts`). Kept out of the hot path so it is swappable.
cluster_seurat <- function(counts, resolution = 1.0, npcs = 30L, nfeatures = 2000L) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("cluster_seurat needs the Seurat package; supply your own `cluster_fn` otherwise.",
         call. = FALSE)
  suppressWarnings(suppressMessages({
    so <- Seurat::CreateSeuratObject(counts = counts)
    so <- Seurat::NormalizeData(so, verbose = FALSE)
    so <- Seurat::FindVariableFeatures(so, nfeatures = nfeatures, verbose = FALSE)
    so <- Seurat::ScaleData(so, verbose = FALSE)
    so <- Seurat::RunPCA(so, npcs = npcs, verbose = FALSE)
    so <- Seurat::FindNeighbors(so, dims = seq_len(npcs), verbose = FALSE)
    so <- Seurat::FindClusters(so, resolution = resolution, verbose = FALSE)
  }))
  paste0("c", as.character(so$seurat_clusters))
}

# Simulate heterotypic doublets between DIFFERENT clusters by summing raw counts.
# `per` is floored at `min_pair` so every included cluster-pair yields enough
# doublets to fit its model (matches run_methodC_reffree*.R). Deterministic given seed.
simulate_cluster_doublets <- function(counts, sing_idx, cl, target_rate, min_pair, seed) {
  set.seed(seed)
  clusters <- sort(unique(cl))
  t2c   <- split(sing_idx, cl)
  n_dbl <- round(length(sing_idx) * target_rate / (1 - target_rate))
  het   <- utils::combn(clusters, 2, simplify = FALSE)
  per   <- max(min_pair, ceiling(n_dbl / length(het)))
  dA <- integer(0); dB <- integer(0); dpl <- character(0)
  for (pr in het) {
    a <- pr[1]; b <- pr[2]
    if (length(t2c[[a]]) < 1 || length(t2c[[b]]) < 1) next
    ca <- sample(t2c[[a]], per, replace = TRUE)
    cb <- sample(t2c[[b]], per, replace = TRUE)
    dA <- c(dA, ca); dB <- c(dB, cb); dpl <- c(dpl, rep(pair_label(a, b), per))
  }
  dcnt <- counts[, dA, drop = FALSE] + counts[, dB, drop = FALSE]
  colnames(dcnt) <- paste0("sim_", seq_len(ncol(dcnt)))
  list(counts = dcnt, pair_label = dpl, clusters = clusters)
}

##############################################################################
# compose_reffree(): reference-free detect + compose from the query alone.
#
# `query`       : gene x cell raw counts.
# `hc_singlets` : logical over cells; round-0 doublet flag is its complement's
#                 removal target. Default: start from all-singlet (rep FALSE) — i.e.
#                 round 0 clusters ALL droplets — matching the scripts.
# `target_rate` : simulated-doublet proportion per round (required; e.g. observed
#                 scDblFinder rate). `cluster_fn`/`score_fn` are injectable.
#
# Returns an `obligato_reffree` object: per-cell detection score, flag (ecdf), and
# composition (top cluster-pair), the final clustering + pair models, and per-round
# diagnostics. Deterministic given `seed`.
##############################################################################
compose_reffree <- function(query, hc_singlets = NULL, target_rate,
                            n_rounds = 4L, resolution = 1.0, npcs = 30L,
                            min_pair = 6L, sing_min = 15L, seed = 1L,
                            cluster_fn = cluster_seurat, score_fn = score_phispace,
                            norm_fn = function(sce) PhiSpace::scranTransf(sce)) {
  stopifnot(is.numeric(target_rate), target_rate > 0, target_rate < 1)
  N <- ncol(query)
  dbl_flag <- if (is.null(hc_singlets)) rep(FALSE, N) else !hc_singlets
  stopifnot(length(dbl_flag) == N)

  rounds <- vector("list", n_rounds + 1L)
  final  <- list()
  for (r in 0:n_rounds) {
    sing_idx <- which(!dbl_flag)
    cl       <- cluster_fn(query[, sing_idx, drop = FALSE], resolution = resolution, npcs = npcs)
    clusters <- sort(unique(cl))

    sim <- simulate_cluster_doublets(query, sing_idx, cl, target_rate, min_pair, seed + r)
    ref <- SingleCellExperiment::SingleCellExperiment(
             assays = list(counts = query[, sing_idx, drop = FALSE]))
    ref <- norm_fn(ref); ref$cluster <- cl

    S_all  <- score_fn(cbind(query, sim$counts), ref, "cluster", sprintf(" reffree-r%d", r))
    S_pool <- S_all[seq_len(N), , drop = FALSE]
    S_dbl  <- S_all[(N + 1L):nrow(S_all), , drop = FALSE]

    Yps <- to_logit(S_pool[sing_idx, , drop = FALSE])
    sing_models <- list()
    for (cc in clusters) { idx <- which(cl == cc)
      if (length(idx) >= sing_min) sing_models[[cc]] <- fit_ln_shrink(Yps[idx, , drop = FALSE]) }
    Ydt <- to_logit(S_dbl); pair_models <- list()
    for (pl in unique(sim$pair_label)) { idx <- which(sim$pair_label == pl)
      if (length(idx) >= min_pair) pair_models[[pl]] <- fit_ln_shrink(Ydt[idx, , drop = FALSE]) }
    if (length(pair_models) == 0L)
      stop(sprintf("No pair models fitted in round %d; lower `resolution` or `min_pair`.", r),
           call. = FALSE)

    F_pool <- compute_features(S_pool, sing_models, pair_models)
    F_dbl  <- compute_features(S_dbl,  sing_models, pair_models)
    bst    <- train_detector(F_pool[sing_idx, , drop = FALSE], F_dbl, seed = seed)
    score_pool <- predict(bst, as.matrix(F_pool))
    thr        <- ecdf_threshold(score_pool, predict(bst, as.matrix(F_dbl)))

    rounds[[r + 1L]] <- list(round = r, n_clean = length(sing_idx),
                             n_clusters = length(clusters), threshold = thr,
                             n_flagged = sum(score_pool > thr))
    final <- list(score = score_pool, threshold = thr, clusters = cl,
                  sing_idx = sing_idx, pair_models = pair_models, S_pool = S_pool)
    dbl_flag <- score_pool > thr
  }

  # Composition: rank every fitted cluster-pair model by likelihood (as the script).
  LL   <- all_pair_ll(final$S_pool, final$pair_models)
  top1 <- max.col(LL, ties.method = "first")
  comp <- data.frame(cell = seq_len(N),
                     top_pair = colnames(LL)[top1],
                     top_loglik = LL[cbind(seq_len(N), top1)],
                     stringsAsFactors = FALSE)

  info <- list(n_cells = N, n_rounds = n_rounds, resolution = resolution, npcs = npcs,
               min_pair = min_pair, target_rate = target_rate, seed = seed,
               threshold = final$threshold, n_final_clusters = length(unique(final$clusters)))
  structure(list(composition = comp, detection_score = final$score,
                 flag = final$score > final$threshold, clusters = final$clusters,
                 pair_models = final$pair_models, rounds = rounds, info = info),
            class = "obligato_reffree")
}
