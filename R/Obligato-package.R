#' Obligato: doublet composition from continuous cell-type scores
#'
#' Obligato identifies the two cell types contributing to a doublet. Its primary
#' interface is [compose()], which can use the built-in detector, external doublet
#' calls, or no detector. [compose_reffree()] provides an iterative clustering mode
#' when a labelled reference is unavailable.
#'
#' @keywords internal
#' @importFrom corpcor cov.shrink
#' @importFrom stats dnorm ecdf predict sd uniroot
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom PhiSpace scranTransf PhiSpace
"_PACKAGE"
