# Write adduct hypotheses as CAMERA-style text

`formatAdducts()` writes the `nBest` most plausible hypotheses of each
feature as text, best first, in CAMERA's format
(`"[M+H]+ 169.037 [M+Na]+ 146.048"`).

## Usage

``` r
formatAdducts(object, ...)

# S4 method for class 'data.frame'
formatAdducts(object, nBest = 3, nFeatures = max(object$feature, 0L), ...)

# S4 method for class 'XcmsExperiment'
formatAdducts(object, nBest = 3, ...)
```

## Arguments

- object:

  A `data.frame` of hypotheses as returned by the `numeric` method of
  [`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md),
  or an annotated
  [xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html).

- ...:

  ignored.

- nBest:

  `numeric(1)` number of hypotheses per feature; `Inf` for all.

- nFeatures:

  `integer(1)`, for the `data.frame` method only: the total number of
  features (length of the result).

## Value

For the `data.frame` method a `character` vector with one entry per
feature (`""` without hypotheses). For the `XcmsExperiment` method the
object with a rewritten `adduct` column; nothing is recomputed.

## Examples

``` r
mz <- c(124.0392, 152.0342, 170.0445, 192.0275, 377.0300)
hints <- annotateAdducts(mz, AdductVoteParam(),
                         intensity = c(81, 50, 100, 11, 2),
                         groups = rep("FG.1", 5))
formatAdducts(hints, nBest = 1, nFeatures = 5)
#> [1] "[M+H-HCOOH]+ 169.037" "[M+H-H2O]+ 169.037"   "[M+H]+ 169.037"      
#> [4] "[M+Na]+ 169.037"      "[2M+K]+ 169.037"     
```
