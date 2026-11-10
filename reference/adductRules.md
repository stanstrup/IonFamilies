# Default adduct and fragment rules

`adductRules` is a list with two `data.frame`s, `positive` and
`negative`, with the adduct and in-source fragment rules from the
'commonMZ' package (`commonMZ::MZ_CAMERA()`), in CAMERA's rule format
and expanded to one to four molecules per ion. The list is long on
purpose: annotations are meant as hints, so a missed true relationship
is worse than a coincidental one. The ranking of
[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
puts the best supported hypothesis first.

## Usage

``` r
data(adductRules)

defaultAdductRules(polarity = c("positive", "negative"))
```

## Format

A list of two `data.frame`s with the columns:

- name:

  rule name, e.g. `"[M+H]+"`.

- nmol:

  number of molecules in the ion.

- charge:

  charge of the ion.

- massdiff:

  mass added to `nmol` neutral molecules to form the ion.

- oidscore:

  CAMERA's rule identifier score.

- quasi:

  1 for quasi-molecular ions, which can anchor a hypothesis.

- ips:

  rule score; higher means more common.

## Source

Generated from 'commonMZ' (MIT licence, by the same author) with
`inst/scripts/make-adductRules.R`.

## Arguments

- polarity:

  `character(1)`, `"positive"` or `"negative"`.

## Value

`defaultAdductRules()`: a `data.frame`.

## Details

`defaultAdductRules()` returns the table for one polarity.

## Examples

``` r
head(defaultAdductRules("positive"))
#>        name nmol charge massdiff oidscore quasi  ips
#> 1    [M+H]+    1      1 1.007276        1     1 1.00
#> 2   [2M+H]+    2      1 1.007276        1     0 1.00
#> 3   [3M+H]+    3      1 1.007276        1     0 1.00
#> 4   [4M+H]+    4      1 1.007276        1     0 1.00
#> 5  [M+2H]2+    1      2 2.014552        2     0 0.75
#> 6 [2M+2H]2+    2      2 2.014552        2     0 0.75
```
