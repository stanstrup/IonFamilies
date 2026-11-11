# Annotate adducts and fragments within feature groups

`annotateAdducts()` finds adducts and in-source fragments by
neutral-mass voting within each feature group (see
[`AdductVoteParam()`](https://stanstrup.github.io/IonFamilies/reference/AdductVoteParam.md)
for the method). Every supported hypothesis is kept, so a feature can
have several ranked alternatives.

## Usage

``` r
annotateAdducts(object, param, ...)

# S4 method for class 'numeric,AdductVoteParam'
annotateAdducts(
  object,
  param,
  intensity,
  groups,
  exclude = rep(FALSE, length(object)),
  ...
)

# S4 method for class 'XcmsExperiment,AdductVoteParam'
annotateAdducts(object, param, ...)
```

## Arguments

- object:

  For the `numeric` method: the m/z values of the features. For the
  `XcmsExperiment` method: an
  [xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html)
  with feature groups.

- param:

  An
  [AdductVoteParam](https://stanstrup.github.io/IonFamilies/reference/AdductVoteParam.md)
  object.

- ...:

  ignored.

- intensity:

  `numeric`, for the `numeric` method only: one representative intensity
  per feature.

- groups:

  for the `numeric` method only: one group label per feature.

- exclude:

  `logical`, for the `numeric` method only: features to leave out,
  normally heavy isotopologues (`isIsotopologue(labels)`). The
  `XcmsExperiment` method uses the `isotopes` column if present.

## Value

For the `numeric` method a `data.frame` with one row per hypothesis and
feature and the columns `feature` (index), `feature_group`, `hypothesis`
(rank within the group), `rank` (rank for this feature, 1 = most
plausible), `M` (neutral mass), `adduct` (rule name), `n_ions` (number
of features supporting `M`), `score` and `intensity`. For the
`XcmsExperiment` method the object, with the columns `adduct` (the best
`nBest` hints as text) and `adduct_hints` (a list column with all
hypotheses per feature) added to
[`featureDefinitions()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html),
and the parameters in `metadata(object)$IonFamilies$adducts`.

## See also

[`AdductVoteParam()`](https://stanstrup.github.io/IonFamilies/reference/AdductVoteParam.md),
[`adductHints()`](https://stanstrup.github.io/IonFamilies/reference/adductHints.md),
[`formatAdducts()`](https://stanstrup.github.io/IonFamilies/reference/formatAdducts.md),
[`bestHypotheses()`](https://stanstrup.github.io/IonFamilies/reference/bestHypotheses.md),
[`annotateIsotopes()`](https://stanstrup.github.io/IonFamilies/reference/annotateIsotopes.md)

## Examples

``` r
mz <- c(95.0124, 124.0392, 152.0342, 170.0445, 171.0496, 192.0275,
        377.0300, 393.0019)
int <- c(71.6, 81.4, 49.8, 100, 4.33, 11.3, 1.76, 1.82)
grp <- rep("FG.1", 8)
iso <- annotateIsotopes(mz, IsotopologueParam(), intensity = int,
                        groups = grp)
hints <- annotateAdducts(mz, AdductVoteParam(), intensity = int,
                         groups = grp, exclude = isIsotopologue(iso))
head(hints)
#>   feature feature_group hypothesis rank         M
#> 1       1          FG.1          2    1 151.02692
#> 2       1          FG.1          3    2 147.05528
#> 3       1          FG.1          4    3  94.00512
#> 4       1          FG.1          9    4 152.01065
#> 5       1          FG.1         10    5 129.04498
#> 6       2          FG.1          1    1 169.03722
#>                                  adduct n_ions score intensity
#> 1 [M+H - O=C=NCH3 (methyl isocyanate)]+      3  2.25    202.80
#> 2                    [M+H-NH3-H2O-H2O]+      3  2.25    182.90
#> 3                                [M+H]+      2  2.00     73.36
#> 4                         [M+H-C2H2O2]+      2  1.50    171.60
#> 5                        [M+H-NH3-H2O]+      2  1.50    121.40
#> 6                          [M+H-HCOOH]+      5  4.25    244.26
```
