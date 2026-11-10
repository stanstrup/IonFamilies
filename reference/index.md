# Package index

## Annotating an XcmsExperiment or plain vectors

Isotope and adduct annotation within feature groups, and access to the
stored results.

- [`annotateIsotopes()`](https://stanstrup.github.io/IonFamilies/reference/annotateIsotopes.md)
  : Annotate isotopes within feature groups
- [`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
  : Annotate adducts and fragments within feature groups
- [`adductHints()`](https://stanstrup.github.io/IonFamilies/reference/adductHints.md)
  : All adduct hypotheses of an annotated object
- [`formatAdducts()`](https://stanstrup.github.io/IonFamilies/reference/formatAdducts.md)
  : Write adduct hypotheses as CAMERA-style text
- [`bestHypotheses()`](https://stanstrup.github.io/IonFamilies/reference/bestHypotheses.md)
  : Reduce adduct hypotheses to one per feature
- [`isIsotopologue()`](https://stanstrup.github.io/IonFamilies/reference/isIsotopologue.md)
  : Is an isotope label a heavy isotopologue?

## Parameters

- [`IsotopologueParam()`](https://stanstrup.github.io/IonFamilies/reference/IsotopologueParam.md)
  [`show(`*`<IsotopologueParam>`*`)`](https://stanstrup.github.io/IonFamilies/reference/IsotopologueParam.md)
  : Parameters for isotope annotation
- [`AdductVoteParam()`](https://stanstrup.github.io/IonFamilies/reference/AdductVoteParam.md)
  [`show(`*`<AdductVoteParam>`*`)`](https://stanstrup.github.io/IonFamilies/reference/AdductVoteParam.md)
  : Parameters for adduct and fragment annotation

## Rules

- [`adductRules`](https://stanstrup.github.io/IonFamilies/reference/adductRules.md)
  [`defaultAdductRules()`](https://stanstrup.github.io/IonFamilies/reference/adductRules.md)
  : Default adduct and fragment rules

## Low-level functions

The building blocks of the adduct annotation, exported for teaching and
testing.

- [`impliedMasses()`](https://stanstrup.github.io/IonFamilies/reference/impliedMasses.md)
  : Neutral masses implied by m/z values and rules
- [`neutralMassVote()`](https://stanstrup.github.io/IonFamilies/reference/neutralMassVote.md)
  : Neutral-mass voting within one feature group

## Package

- [`IonFamilies`](https://stanstrup.github.io/IonFamilies/reference/IonFamilies-package.md)
  [`IonFamilies-package`](https://stanstrup.github.io/IonFamilies/reference/IonFamilies-package.md)
  : IonFamilies: isotope and adduct annotation of LC-MS feature groups
