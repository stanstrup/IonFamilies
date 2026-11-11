# Reduce adduct hypotheses to one per feature

`bestHypotheses()` reproduces CAMERA's final choice: within each group
it accepts hypotheses from the best ranked down, skipping any hypothesis
that contains a feature already explained by an accepted one. Each
feature then keeps at most one neutral mass.

## Usage

``` r
bestHypotheses(x)
```

## Arguments

- x:

  `data.frame` of hypotheses as returned by
  [`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
  or
  [`adductHints()`](https://stanstrup.github.io/IonFamilies/reference/adductHints.md).

## Value

The accepted rows of `x`.

## Examples

``` r
mz <- c(124.0392, 152.0342, 170.0445, 192.0275, 377.0300)
hints <- annotateAdducts(mz, AdductVoteParam(),
                         intensity = c(81, 50, 100, 11, 2),
                         groups = rep("FG.1", 5))
bestHypotheses(hints)
#>   feature feature_group hypothesis rank        M       adduct n_ions score
#> 1       1          FG.1          1    1 169.0372 [M+H-HCOOH]+      5  4.25
#> 2       2          FG.1          1    1 169.0372   [M+H-H2O]+      5  4.25
#> 3       3          FG.1          1    1 169.0372       [M+H]+      5  4.25
#> 4       4          FG.1          1    1 169.0372      [M+Na]+      5  4.25
#> 5       5          FG.1          1    1 169.0372      [2M+K]+      5  4.25
#>   intensity
#> 1       244
#> 2       244
#> 3       244
#> 4       244
#> 5       244
```
