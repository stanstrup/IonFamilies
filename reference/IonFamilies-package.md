# IonFamilies: isotope and adduct annotation of LC-MS feature groups

One compound gives several LC-MS features: its quasi-molecular ion,
other adducts, in-source fragments, and the isotopes of all of these.
This package annotates that *ion family* within groups of features,
typically the groups made with
[`MsFeatures::groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
on an
[xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html).

## Details

The workflow has two steps:

1.  [`annotateIsotopes()`](https://stanstrup.github.io/IonFamilies/reference/annotateIsotopes.md)
    with an
    [`IsotopologueParam()`](https://stanstrup.github.io/IonFamilies/reference/IsotopologueParam.md)
    finds isotope patterns with
    [`MetaboCoreUtils::isotopologues()`](https://rdrr.io/pkg/MetaboCoreUtils/man/isotopologues.html)
    and labels them in CAMERA syntax (`[3][M]+`, `[3][M+1]+`).

2.  [`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
    with an
    [`AdductVoteParam()`](https://stanstrup.github.io/IonFamilies/reference/AdductVoteParam.md)
    finds adducts and fragments by neutral-mass voting with a rule list
    ([`defaultAdductRules()`](https://stanstrup.github.io/IonFamilies/reference/adductRules.md)).
    All supported hypotheses are kept and ranked
    ([`adductHints()`](https://stanstrup.github.io/IonFamilies/reference/adductHints.md));
    the best `nBest` are written as text
    ([`formatAdducts()`](https://stanstrup.github.io/IonFamilies/reference/formatAdducts.md)).

Annotations are *hints*: they suggest what to check, they do not
identify a compound.

## See also

Useful links:

- <https://stanstrup.github.io/IonFamilies>

- <https://github.com/stanstrup/IonFamilies>

- Report bugs at <https://github.com/stanstrup/IonFamilies/issues>

## Author

**Maintainer**: Jan Stanstrup <stanstrup@gmail.com>

Authors:

- Jan Stanstrup <stanstrup@gmail.com>
