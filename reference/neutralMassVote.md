# Neutral-mass voting within one feature group

`neutralMassVote()` is the core of
[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
for a single group of features. Every feature is tried as every
quasi-molecular ion (`quasi == 1`). A neutral mass is a hypothesis when
at least one other feature is explained by some rule as an ion of the
same mass. Duplicate hypotheses (same features with the same rules) are
merged, and the rest are ranked by summed rule score (`ips`), then by
explained intensity.

## Usage

``` r
neutralMassVote(mz, intensity, rules, ppm = 10, mzabs = 0.005)
```

## Arguments

- mz:

  `numeric` m/z values of the group's features.

- intensity:

  `numeric` representative intensities of the features.

- rules:

  `data.frame` with rules in CAMERA's format.

- ppm:

  `numeric(1)` relative tolerance in ppm.

- mzabs:

  `numeric(1)` absolute tolerance.

## Value

`data.frame` with one row per hypothesis and member feature and the
columns `feature` (index into `mz`), `hypothesis` (rank), `M`, `adduct`,
`n_ions`, `score` and `intensity`.

## Examples

``` r
mz <- c(124.0392, 152.0342, 170.0445, 192.0275, 377.0300)
neutralMassVote(mz, c(81, 50, 100, 11, 2), defaultAdductRules())
#>    feature hypothesis        M           adduct n_ions score intensity
#> 1        3          1 169.0372           [M+H]+      5  4.25       244
#> 2        4          1 169.0372          [M+Na]+      5  4.25       244
#> 3        5          1 169.0372          [2M+K]+      5  4.25       244
#> 4        1          1 169.0372     [M+H-HCOOH]+      5  4.25       244
#> 5        2          1 169.0372       [M+H-H2O]+      5  4.25       244
#> 6        3          2 187.0478       [M+H-H2O]+      3  1.75       231
#> 7        2          2 187.0478    [M+H-(H2O)2]+      3  1.75       231
#> 8        1          2 187.0478 [M+H-H2O-HCOOH]+      3  1.75       231
#> 9        2          3 151.0269           [M+H]+      2  1.75       131
#> 10       1          3 151.0269        [M+H-CO]+      2  1.75       131
#> 11       3          4 147.0553          [M+Na]+      2  1.75       111
#> 12       4          4 147.0553       [M+2Na-H]+      2  1.75       111
#> 13       5          5 338.0668           [M+K]+      2  1.75       102
#> 14       3          5 338.0668         [M+2H]2+      2  1.75       102
#> 15       1          6 101.0500          [M+Na]+      2  1.50       181
#> 16       3          6 101.0500    [M+H+CHOONa]+      2  1.50       181
#> 17       3          7 131.0813           [M+K]+      2  1.50       111
#> 18       4          7 131.0813      [M+Na+K-H]+      2  1.50       111
#> 19       1          8 123.0319           [M+H]+      2  1.50        92
#> 20       4          8 123.0319    [M+H+CHOONa]+      2  1.50        92
#> 21       4          9 209.0308       [M+H-H2O]+      2  1.25        61
#> 22       2          9 209.0308    [M+H-C2H2O2]+      2  1.25        61
```
