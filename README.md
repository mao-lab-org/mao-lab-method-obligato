# Obligato

**Doublet composition identification from continuous cell-type scores.**

`Obligato` identifies the cell-type *composition* of doublets in single-cell RNA-seq
data — which two cell types a doublet contains, not merely that a droplet is abnormal.
It fits per-pair logit-normal models over [PhiSpace](https://github.com/jiadongm/PhiSpace)
continuous cell-type scores, treating within-lineage co-elevation as a baseline rather than
as signal. Composition can be applied downstream of any doublet detector; a detector is also
provided.

The name is the anglicised variant of the musical *obbligato* — a part that cannot be
omitted. Detection tells you a droplet is a doublet; composition is the part you cannot leave
out if the removal decision is to be informed.

> Early development (`0.0.0.9000`). API may change before the first tagged release.

## Installation

`Obligato` depends on PhiSpace, which is distributed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("jiadongm/PhiSpace/pkg")   # PhiSpace
remotes::install_github("jiadongm/Obligato")       # this package
```

## Usage

`compose()` is the main entry point. Detection is an optional module.

```r
library(Obligato)

## Full pipeline: our detector + composition
res <- compose(query, reference, phenotypes = "celltype",
               hc_singlets = hc_flag)          # detector = "builtin" (default)
res$composition   # per-cell c1 / c2
res$flag          # doublet calls
res$detection_score

## Bring your own detector — composition only
res <- compose(query, reference, phenotypes = "celltype",
               hc_singlets = hc_flag, detector = my_scDblFinder_calls)

## No detection, compose everything
res <- compose(query, reference, phenotypes = "celltype",
               hc_singlets = hc_flag, detector = "none")
```

## Evaluation

Exported metrics for benchmarking composition against ground truth, including the
maximum achievable accuracy of a cluster-output grouping:

```r
composition_accuracy(pred_A, pred_B, true_A, true_B, lineage = lin)
max_achievable_accuracy(true_A, true_B, decode, definition = "vocabulary")
```

## License

MIT © Jiadong Mao, George Howitt, Michelle Meier
