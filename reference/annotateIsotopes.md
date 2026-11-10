# Annotate isotopes within feature groups

`annotateIsotopes()` labels isotope patterns within each feature group
in CAMERA syntax: `"[3][M]+"` is the monoisotopic peak of pattern 3,
`"[3][M+1]+"` its first isotopologue, and so on. Features that are not
part of a pattern get `""`.

## Usage

``` r
annotateIsotopes(object, param, ...)

# S4 method for class 'numeric,IsotopologueParam'
annotateIsotopes(object, param, intensity, groups, ...)

# S4 method for class 'XcmsExperiment,IsotopologueParam'
annotateIsotopes(object, param, ...)
```

## Arguments

- object:

  For the `numeric` method: the m/z values of the features. For the
  `XcmsExperiment` method: an
  [xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html)
  with feature groups (see
  [`MsFeatures::groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)).

- param:

  An
  [IsotopologueParam](https://stanstrup.github.io/IonFamilies/reference/IsotopologueParam.md)
  object.

- ...:

  ignored.

- intensity:

  `numeric`, for the `numeric` method only: one representative intensity
  per feature, e.g. the median over the study samples.

- groups:

  for the `numeric` method only: one group label per feature.

## Value

For the `numeric` method a `character` vector with one label per
feature. For the `XcmsExperiment` method the object, with the labels in
the `isotopes` column of
[`featureDefinitions()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
and the parameters in `metadata(object)$IonFamilies$isotopes`.

## See also

[`IsotopologueParam()`](https://stanstrup.github.io/IonFamilies/reference/IsotopologueParam.md),
[`isIsotopologue()`](https://stanstrup.github.io/IonFamilies/reference/isIsotopologue.md),
[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)

## Examples

``` r
## Ions of 2-furoylglycine (M = 169.0375) in one feature group
mz <- c(95.0124, 124.0392, 152.0342, 170.0445, 171.0496, 192.0275,
        377.0300, 393.0019)
int <- c(71.6, 81.4, 49.8, 100, 4.33, 11.3, 1.76, 1.82)
annotateIsotopes(mz, IsotopologueParam(), intensity = int,
                 groups = rep("FG.1", 8))
#> [1] ""          ""          ""          "[1][M]+"   "[1][M+1]+" ""         
#> [7] ""          ""         
```
