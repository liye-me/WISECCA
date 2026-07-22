
# WISECCA

![WISECCA workflow](man/figures/WISECCA_workflow.png)

WISECCA implements **Weighted co-expression network-Informed Sparse
intEgrative Canonical Correlation Analysis** for integrating molecular
omics data with neuroimaging features.

The workflow combines weighted co-expression network analysis, sparse
canonical correlation analysis, animal-level cross-validation,
block-bootstrap stability analysis, LASSO regression, and mixed-effects
modeling.

## Installation

Install the development version from GitHub:

``` r
remotes::install_github("liye-me/WISECCA")
```

Load the package:

``` r
library(WISECCA)
```

## Input data

WISECCA requires three matched objects:

- `omics`: an observation-by-feature molecular matrix;
- `imaging`: an observation-by-feature imaging matrix;
- `meta`: a sample metadata data frame.

The metadata should contain the following columns:

- `sample_id`
- `animal`
- `brain`
- `age`

The rows of `omics`, `imaging`, and `meta` must represent the same
observations in the same order.

## Main workflow

``` r
result <- wisecca(
  omics = omics,
  imaging = imaging,
  meta = meta,
  id_col = "sample_id",
  animal_col = "animal",
  brain_col = "brain",
  age_col = "age",
  power = 6,
  run_bootstrap = TRUE,
  bootstrap_args = list(
    n_boot = 500
  )
)
```

The main workflow includes:

1.  Input validation
2.  WGCNA module construction
3.  Animal-level sparse CCA cross-validation
4.  Final sparse CCA fitting
5.  Animal-block bootstrap analysis

## Main functions

| Function | Purpose |
|----|----|
| `wise_check_input()` | Validate and match input data |
| `wise_modules()` | Construct weighted co-expression modules |
| `wise_scca()` | Fit a sparse CCA model |
| `wise_cv_scca()` | Select sparse CCA penalties using animal-level cross-validation |
| `wise_bootstrap()` | Evaluate feature-selection stability |
| `wise_lasso()` | Fit animal-level cross-validated LASSO models |
| `wise_mixed()` | Fit mixed-effects models |
| `wisecca()` | Run the complete workflow |

## Visualization

``` r
plot(result, type = "scores")
plot(result, type = "cv")
plot(result, type = "weights", side = "x")
plot(result, type = "bootstrap", side = "z")
```

Optional results can be plotted after enabling the corresponding
analyses:

``` r
plot(result, type = "lasso")
plot(result, type = "mixed")
```

## Development status

WISECCA is under active development. The package structure and core
analysis functions are available, while complete end-to-end validation
on simulated and external datasets is still in progress.

## License

MIT License.
