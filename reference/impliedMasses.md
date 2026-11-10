# Neutral masses implied by m/z values and rules

For an ion with m/z `mz` formed by a rule with `nmol` molecules, charge
`charge` and mass shift `massdiff`, the neutral mass is
`(abs(charge) * mz - massdiff) / nmol`. `impliedMasses()` computes this
for every combination of feature and rule. It is the basic operation of
rule-based adduct annotation and is exported for teaching and testing.

## Usage

``` r
impliedMasses(mz, rules)
```

## Arguments

- mz:

  `numeric` m/z values.

- rules:

  `data.frame` with rules in CAMERA's format (see
  [adductRules](https://stanstrup.github.io/IonFamilies/reference/adductRules.md)).

## Value

`matrix` with one row per m/z value and one column per rule.

## Examples

``` r
rules <- defaultAdductRules("positive")[1:3, ]
impliedMasses(c(170.0445, 192.0275), rules)
#>        [M+H]+  [2M+H]+  [3M+H]+
#> [1,] 169.0372 84.51861 56.34574
#> [2,] 191.0202 95.51011 63.67341
```
