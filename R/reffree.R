# Reference-free composition iterates detection, removal, and reclustering.
# Cluster identities act as provisional types, allowing the same likelihood model
# to operate without an external labelled reference.

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

# Simulate heterotypic doublets between clusters. The per-pair count is floored
# at `min_pair` so every included pair can be fitted.
simulate_cluster_doublets <- function(counts, sing_idx, cl, target_rate,
                                      min_pair, seed) {
  if (length(sing_idx) != length(cl) || anyNA(sing_idx) || anyNA(cl)) {
    stop("Cluster labels must contain one non-missing value per singlet.",
         call. = FALSE)
  }
  clusters <- sort(unique(cl))
  if (length(clusters) < 2L) {
    stop("Reference-free composition requires at least two clusters.", call. = FALSE)
  }
  min_pair <- .validate_scalar_integer(min_pair, "min_pair", 1L)
  .with_seed(seed, {
    t2c <- split(sing_idx, cl)
    n_dbl <- round(length(sing_idx) * target_rate / (1 - target_rate))
    het <- utils::combn(clusters, 2L, simplify = FALSE)
    per <- max(min_pair, ceiling(n_dbl / length(het)))
    dA <- integer(); dB <- integer(); dpl <- character()
    for (pr in het) {
      a <- pr[1L]; b <- pr[2L]
      ca <- sample(t2c[[a]], per, replace = TRUE)
      cb <- sample(t2c[[b]], per, replace = TRUE)
      dA <- c(dA, ca); dB <- c(dB, cb)
      dpl <- c(dpl, rep(pair_label(a, b), per))
    }
    dcnt <- counts[, dA, drop = FALSE] + counts[, dB, drop = FALSE]
    colnames(dcnt) <- paste0("sim_", seq_len(ncol(dcnt)))
    list(counts = dcnt, pair_label = dpl, clusters = clusters)
  })
}

#' Reference-free iterative doublet detection and composition
#'
#' Repeatedly clusters the currently unflagged cells, simulates cross-cluster
#' doublets, detects doublets, and removes them before the next clustering round.
#' Cluster pairs are proxies for cell-type pairs. This mode currently uses
#' PhiSpace score features only; the reference-based library-size detector is not
#' silently substituted because it has not been validated for iterative clusters.
#'
#' @param query A genes-by-cells matrix-like object of raw non-negative counts.
#' @param hc_singlets Optional initial non-missing logical singlet calls.
#' @param target_rate Expected doublet fraction, strictly between zero and one.
#' @param n_rounds Number of removal/reclustering updates after round zero.
#' @param resolution,npcs Clustering parameters passed to `cluster_fn`.
#' @param min_pair Minimum simulated cells required for a cluster-pair model.
#' @param sing_min Minimum cells required for a cluster singlet model.
#' @param seed Non-negative integer controlling clustering-time scoring, simulation,
#'   and detector training.
#' @param cluster_fn Function returning one cluster label per supplied cell.
#' @param score_fn Advanced scoring function with the same interface as
#'   `score_phispace()`.
#' @param norm_fn Function that normalizes the temporary reference object.
#' @return An `obligato_reffree` object with per-cell composition, scores,
#'   flags, final cluster membership, round diagnostics, and provenance.
#' @export
compose_reffree <- function(query, hc_singlets = NULL, target_rate,
                            n_rounds = 4L, resolution = 1.0, npcs = 30L,
                            min_pair = 6L, sing_min = 15L, seed = 1L,
                            cluster_fn = cluster_seurat, score_fn = score_phispace,
                            norm_fn = function(sce) PhiSpace::scranTransf(sce)) {
  .validate_counts(query)
  if (length(target_rate) != 1L || !is.numeric(target_rate) ||
      is.na(target_rate) || target_rate <= 0 || target_rate >= 1) {
    stop("`target_rate` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  n_rounds <- .validate_scalar_integer(n_rounds, "n_rounds", 0L)
  npcs <- .validate_scalar_integer(npcs, "npcs", 1L)
  min_pair <- .validate_scalar_integer(min_pair, "min_pair", 2L)
  sing_min <- .validate_scalar_integer(sing_min, "sing_min", 2L)
  seed <- .validate_scalar_integer(seed, "seed", 0L)
  if (length(resolution) != 1L || !is.numeric(resolution) ||
      is.na(resolution) || resolution <= 0) {
    stop("`resolution` must be a single positive number.", call. = FALSE)
  }

  N <- ncol(query)
  if (is.null(hc_singlets)) {
    dbl_flag <- rep(FALSE, N)
  } else {
    if (!is.logical(hc_singlets) || length(hc_singlets) != N ||
        anyNA(hc_singlets)) {
      stop("`hc_singlets` must contain one non-missing logical value per cell.",
           call. = FALSE)
    }
    dbl_flag <- !hc_singlets
  }

  rounds <- vector("list", n_rounds + 1L)
  final <- list()
  for (r in 0:n_rounds) {
    sing_idx <- which(!dbl_flag)
    if (length(sing_idx) < 2L) {
      stop(sprintf("Fewer than two unflagged cells remain before round %d.", r),
           call. = FALSE)
    }
    cl <- .with_seed(seed + r,
      cluster_fn(query[, sing_idx, drop = FALSE],
                 resolution = resolution, npcs = npcs))
    if (length(cl) != length(sing_idx) || anyNA(cl)) {
      stop(sprintf("The clustering function returned invalid labels in round %d.", r),
           call. = FALSE)
    }
    cl <- as.character(cl)
    clusters <- sort(unique(cl))
    if (length(clusters) < 2L) {
      stop(sprintf("Reference-free mode found fewer than two clusters in round %d.", r),
           call. = FALSE)
    }

    sim <- simulate_cluster_doublets(
      query, sing_idx, cl, target_rate, min_pair, seed + r)
    ref <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = query[, sing_idx, drop = FALSE]))
    ref <- .with_seed(seed + 50L + r, norm_fn(ref))
    ref$cluster <- cl

    score_counts <- cbind(query, sim$counts)
    S_all <- .with_seed(seed + 100L + r,
      score_fn(score_counts, ref, "cluster", sprintf(" reffree-r%d", r)))
    S_all <- .validate_scores(
      S_all, ncol(score_counts), sprintf("Reference-free scores in round %d", r))
    S_pool <- S_all[seq_len(N), , drop = FALSE]
    S_dbl <- S_all[N + seq_len(ncol(sim$counts)), , drop = FALSE]

    Yps <- to_logit(S_pool[sing_idx, , drop = FALSE])
    sing_models <- list()
    for (cc in clusters) {
      idx <- which(cl == cc)
      if (length(idx) >= sing_min) {
        sing_models[[cc]] <- fit_ln_shrink(Yps[idx, , drop = FALSE])
      }
    }
    Ydt <- to_logit(S_dbl)
    pair_models <- list()
    for (pl in unique(sim$pair_label)) {
      idx <- which(sim$pair_label == pl)
      if (length(idx) >= min_pair) {
        pair_models[[pl]] <- fit_ln_shrink(Ydt[idx, , drop = FALSE])
      }
    }
    if (!length(sing_models) || !length(pair_models)) {
      stop(sprintf(
        "No usable singlet or pair models in round %d; lower clustering resolution or model minima.",
        r), call. = FALSE)
    }

    F_pool <- compute_features(S_pool, sing_models, pair_models)
    F_dbl <- compute_features(S_dbl, sing_models, pair_models)
    bst <- train_detector(
      F_pool[sing_idx, , drop = FALSE], F_dbl, seed = seed + r)
    score_pool <- as.numeric(predict(bst, as.matrix(F_pool)))
    sim_score <- as.numeric(predict(bst, as.matrix(F_dbl)))
    threshold <- ecdf_threshold(score_pool, sim_score)

    rounds[[r + 1L]] <- list(
      round = r, n_clean = length(sing_idx), n_clusters = length(clusters),
      threshold = threshold, n_flagged = sum(score_pool > threshold))
    cluster_map <- rep(NA_character_, N)
    cluster_map[sing_idx] <- cl
    final <- list(
      score = score_pool, threshold = threshold, clusters = cluster_map,
      sing_idx = sing_idx, pair_models = pair_models, S_pool = S_pool,
      classifier = bst, feature_names = names(F_pool))
    dbl_flag <- score_pool > threshold
  }

  LL <- all_pair_ll(final$S_pool, final$pair_models)
  top1 <- max.col(LL, ties.method = "first")
  comp <- data.frame(
    cell = seq_len(N), top_pair = colnames(LL)[top1],
    top_loglik = LL[cbind(seq_len(N), top1)], stringsAsFactors = FALSE)

  info <- list(
    n_cells = N, n_rounds = n_rounds, resolution = resolution, npcs = npcs,
    min_pair = min_pair, sing_min = sing_min, target_rate = target_rate,
    seed = seed, threshold = final$threshold,
    detection_features = "PhiSpace scores",
    n_final_clusters = length(unique(final$clusters[!is.na(final$clusters)])))
  structure(
    list(composition = comp, detection_score = final$score,
         flag = final$score > final$threshold, clusters = final$clusters,
         pair_models = final$pair_models,
         detection_model = list(
           classifier = final$classifier, threshold = final$threshold,
           feature_names = final$feature_names),
         rounds = rounds, info = info, provenance = .package_provenance()),
    class = "obligato_reffree")
}
