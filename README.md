# Obligato

**Doublet composition identification from continuous cell-type scores.**

`Obligato` identifies which two cell types contribute to a doublet in single-cell
RNA-seq data. It fits per-pair logit-normal models over
[PhiSpace](https://github.com/jiadongm/PhiSpace) scores. Composition can be used
with Obligato's built-in detector, calls from another detector, or no detector.

The name is the anglicised variant of the musical *obbligato*: a part that cannot
be omitted. Detection says that a droplet is likely to be a doublet; composition
provides the information needed to decide what to do with it.

> Early development (`0.0.0.9001`). The API may change before the first tagged release.

## Installation

```r
# install.packages("BiocManager")
BiocManager::install(c("remotes", "mao-lab-org/mao-lab-method-obligato"))
```

This installs Obligato together with its CRAN, Bioconductor, and GitHub
dependencies, including PhiSpace.

## Inputs

The primary function, `compose()`, expects:

- `query`: a genes-by-cells raw count matrix;
- `reference`: a labelled `SingleCellExperiment` with a `logcounts` assay;
- `phenotypes`: the reference column containing cell-type labels; and
- `hc_singlets`: one logical value per query cell identifying trusted singlets
  used to learn the models.

Optional `hc_labels` may contain one label per query cell or one label per
high-confidence singlet. If omitted, labels are obtained from each singlet's
highest PhiSpace score.

## Reference-based analysis

```r
library(Obligato)

res <- compose(
  query,
  reference,
  phenotypes = "celltype",
  hc_singlets = hc_flag
)

res$composition       # one candidate pair per query cell
res$detection_score   # built-in detector score
res$flag              # calls at the fitted ECDF crossover threshold

# Analyse compositions among called doublets.
called_composition <- res$composition[res$flag[res$composition$cell], ]
```

The built-in detector combines PhiSpace score features with total counts,
detected genes, and a cell-type-relative library-size feature. The composition
model itself remains based on PhiSpace scores.

### Bring your own detector

Logical calls are accepted directly:

```r
res <- compose(
  query, reference, "celltype", hc_flag,
  detector = my_doublet_calls
)
```

Numeric detector scores require an explicit threshold:

```r
res <- compose(
  query, reference, "celltype", hc_flag,
  detector = my_doublet_scores,
  detector_threshold = 0.35
)
```

Use `detector = "none"` to perform composition without making doublet calls.

### Multiple annotation levels

The first annotation level defines the reported composition. Additional levels
contribute detection features:

```r
res <- compose(
  query,
  reference = list(reference_l1, reference_l2),
  phenotypes = c("celltype_l1", "celltype_l2"),
  hc_singlets = hc_flag
)
```

### Alternative composition outputs

```r
top3 <- compose(
  query, reference, "celltype", hc_flag,
  detector = "none", output = "topk", k = 3
)

sets <- compose(
  query, reference, "celltype", hc_flag,
  detector = "none", output = "conformal", alpha = 0.10
)
```

Conformal sets are calibrated using independent synthetic doublets. Their nominal
coverage relies on calibration and target doublets being exchangeable; users
should verify empirical coverage when ground truth is available.

## Reference-free analysis

`compose_reffree()` iteratively clusters unflagged cells and reports cluster-pair
composition when no labelled reference is available:

```r
res_rf <- compose_reffree(query, target_rate = 0.10)
```

This mode currently uses PhiSpace score features for detection. It does not yet
use the library-size feature block validated for the reference-based detector.

## Evaluation

```r
acc <- composition_accuracy(
  pred_A, pred_B, true_A, true_B, lineage = lineage
)

decode <- majority_decode(singlet_group, singlet_type)
ceiling <- max_achievable_accuracy(
  true_A, true_B, decode, definition = "vocabulary"
)

pr <- compute_auprc(is_doublet, detection_score)
```

Missing composition predictions count as incorrect. `composition_accuracy()`
returns both overall accuracy and accuracy balanced equally across observed true
pairs.

Use `help(package = "Obligato")` or `?compose` for the installed documentation.

## License

GNU Affero General Public License v3 (AGPL-3) © Jiadong Mao
