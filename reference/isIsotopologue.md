# Is an isotope label a heavy isotopologue?

`isIsotopologue()` returns `TRUE` for labels of heavy isotopologues
(`"[3][M+1]+"`, `"[3][M+2]-"`) and `FALSE` for monoisotopic peaks
(`"[3][M]+"`), adduct labels such as `"[M+H]+"`, and `""`.

## Usage

``` r
isIsotopologue(x)
```

## Arguments

- x:

  `character` with isotope labels as written by
  [`annotateIsotopes()`](https://stanstrup.github.io/IonFamilies/reference/annotateIsotopes.md)
  or CAMERA.

## Value

`logical` of the same length as `x`.

## Examples

``` r
isIsotopologue(c("[1][M]+", "[1][M+1]+", "[M+H]+", ""))
#> [1] FALSE  TRUE FALSE FALSE
```
