# All adduct hypotheses of an annotated object

`adductHints()` returns all neutral-mass hypotheses stored by
[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
as one long `data.frame` (one row per hypothesis and feature, columns as
described in
[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)).

## Usage

``` r
adductHints(object, ...)

# S4 method for class 'XcmsExperiment'
adductHints(object, ...)
```

## Arguments

- object:

  An
  [xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html)
  annotated with
  [`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md).

- ...:

  ignored.

## Value

A `data.frame`.

## Examples

``` r
library(xcms)
#> Loading required package: BiocParallel
#> 
#> This is xcms version 4.10.1 
library(MsFeatures)
if (requireNamespace("faahKO", quietly = TRUE)) {
    x <- loadXcmsData("xmse")
    x <- groupFeatures(x, SimilarRtimeParam(diffRt = 10))
    x <- annotateAdducts(x, AdductVoteParam())
    head(adductHints(x))
}
#> No isotope annotation found; consider running 'annotateIsotopes()' first so that heavy isotopologues are not used as adducts.
#>   feature feature_group hypothesis rank        M              adduct n_ions
#> 1       2        FG.005          2    1 263.0662        [M+H-C3H9N]+      2
#> 2       2        FG.005          3    2 296.0368 [M+H-C2H4O2-CH3OH]+      2
#> 3       2        FG.005          5    3 296.0033  [M+H-HCOOH-HCOOH]+      2
#> 4       3        FG.005          1    1 240.0368      [M+H-NH3-H2O]+      2
#> 5       3        FG.005          4    2 275.0368        [M+H-C4H6O]+      2
#> 6       4        FG.078          1    1 301.1033  [M+H-S-NH3-HCOOH]+      2
#>   score  intensity
#> 1  1.50  1948834.6
#> 2  1.50  1932821.3
#> 3  1.25 14730092.5
#> 4  1.50 13328259.5
#> 5  1.50   501911.2
#> 6  1.25  1093366.1
```
