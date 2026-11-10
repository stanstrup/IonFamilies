# Parameters for adduct and fragment annotation

`AdductVoteParam()` creates the parameter object for
[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md).

## Usage

``` r
AdductVoteParam(
  rules = NULL,
  ppm = 10,
  mzabs = 0.005,
  polarity = c("positive", "negative"),
  subset = integer(),
  value = "into",
  nBest = 3
)

# S4 method for class 'AdductVoteParam'
show(object)
```

## Arguments

- rules:

  `data.frame` with the adduct and fragment rules in CAMERA's format,
  with columns `name`, `nmol`, `charge`, `massdiff`, `quasi` and `ips`.
  `NULL` (the default) uses
  [`defaultAdductRules()`](https://stanstrup.github.io/IonFamilies/reference/adductRules.md)
  for `polarity`.

- ppm:

  `numeric(1)` relative tolerance (in ppm) for two neutral masses to be
  considered equal.

- mzabs:

  `numeric(1)` absolute tolerance added to `ppm`.

- polarity:

  `character(1)`, `"positive"` or `"negative"`. Selects the default
  rules.

- subset:

  `integer` with the indices of the samples used to compute the
  representative intensity of each feature. An empty vector (the
  default) uses all samples. Only used by the `XcmsExperiment` method.

- value:

  `character(1)` the type of feature value used for the intensities, see
  [`xcms::featureValues()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-peak-grouping-results.html).
  Only used by the `XcmsExperiment` method.

- nBest:

  `numeric(1)` how many hypotheses per feature are written to the
  `adduct` text column. `Inf` writes all of them. All hypotheses are
  stored regardless; see
  [`adductHints()`](https://stanstrup.github.io/IonFamilies/reference/adductHints.md)
  and
  [`formatAdducts()`](https://stanstrup.github.io/IonFamilies/reference/formatAdducts.md).

- object:

  an `AdductVoteParam` object.

## Value

An `AdductVoteParam` object.

## Details

Within each feature group, every feature is tried as every
quasi-molecular ion of the rule list (rules with `quasi == 1`, such as
`[M+H]+`). The neutral mass this implies is a hypothesis if at least one
other feature of the group is explained by some rule as an ion of the
same neutral mass. All hypotheses are kept and ranked by the summed rule
score (`ips`) and then by the summed intensity of the features they
explain.

## See also

[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md),
[`defaultAdductRules()`](https://stanstrup.github.io/IonFamilies/reference/adductRules.md),
[`IsotopologueParam()`](https://stanstrup.github.io/IonFamilies/reference/IsotopologueParam.md)

## Examples

``` r
AdductVoteParam()
#> Object of class "AdductVoteParam"
#>  rules: 159 rules, 5 quasi-molecular
#>  ppm: 10 
#>  mzabs: 0.005 
#>  polarity: positive 
#>  value: into 
#>  subset: all samples 
#>  nBest: 3 
AdductVoteParam(polarity = "negative", nBest = Inf)
#> Object of class "AdductVoteParam"
#>  rules: 147 rules, 5 quasi-molecular
#>  ppm: 10 
#>  mzabs: 0.005 
#>  polarity: negative 
#>  value: into 
#>  subset: all samples 
#>  nBest: Inf 
```
