# Parameters for isotope annotation

`IsotopologueParam()` creates the parameter object for
[`annotateIsotopes()`](https://stanstrup.github.io/IonFamilies/reference/annotateIsotopes.md).
Isotope patterns are searched within each feature group with
[`MetaboCoreUtils::isotopologues()`](https://rdrr.io/pkg/MetaboCoreUtils/man/isotopologues.html),
which checks both the mass differences and the intensity ratios of the
candidate peaks.

## Usage

``` r
IsotopologueParam(
  ppm = 10,
  tolerance = 0.005,
  polarity = c("positive", "negative"),
  subset = integer(),
  value = "into"
)

# S4 method for class 'IsotopologueParam'
show(object)
```

## Arguments

- ppm:

  `numeric(1)` relative m/z tolerance (in ppm) passed to
  [`MetaboCoreUtils::isotopologues()`](https://rdrr.io/pkg/MetaboCoreUtils/man/isotopologues.html).

- tolerance:

  `numeric(1)` absolute m/z tolerance passed to
  [`MetaboCoreUtils::isotopologues()`](https://rdrr.io/pkg/MetaboCoreUtils/man/isotopologues.html).

- polarity:

  `character(1)`, `"positive"` or `"negative"`. Only used for the
  `+`/`-` suffix of the labels.

- subset:

  `integer` with the indices of the samples used to compute the
  representative intensity of each feature (the median over these
  samples). An empty vector (the default) uses all samples. Only used by
  the `XcmsExperiment` method.

- value:

  `character(1)` the type of feature value used for the intensities, see
  [`xcms::featureValues()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-peak-grouping-results.html).
  Only used by the `XcmsExperiment` method.

- object:

  an `IsotopologueParam` object.

## Value

An `IsotopologueParam` object.

## See also

[`annotateIsotopes()`](https://stanstrup.github.io/IonFamilies/reference/annotateIsotopes.md),
[`AdductVoteParam()`](https://stanstrup.github.io/IonFamilies/reference/AdductVoteParam.md)

## Examples

``` r
IsotopologueParam()
#> Object of class "IsotopologueParam"
#>  ppm: 10 
#>  tolerance: 0.005 
#>  polarity: positive 
#>  value: into 
#>  subset: all samples 
IsotopologueParam(ppm = 5, polarity = "negative")
#> Object of class "IsotopologueParam"
#>  ppm: 5 
#>  tolerance: 0.005 
#>  polarity: negative 
#>  value: into 
#>  subset: all samples 
```
