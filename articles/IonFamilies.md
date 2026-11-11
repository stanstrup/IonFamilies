# Getting started with IonFamilies

## Introduction

One compound in an LC-MS experiment gives many features: its protonated
or deprotonated molecule, other adducts such as `[M+Na]+`, in-source
fragments such as `[M+H-H2O]+`, and the isotopes of all of these.
Together they form the compound’s *ion family*. Recognising ion families
reduces redundancy and helps to work out the neutral mass of unknown
compounds.

[`groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
from the *MsFeatures* and *xcms* packages groups features that probably
come from the same compound (similar retention time, correlated
abundances across samples, similar peak shapes). It deliberately stops
there. *IonFamilies* takes the groups and annotates, within each group,

1.  **isotope patterns**, with
    [`MetaboCoreUtils::isotopologues()`](https://rdrr.io/pkg/MetaboCoreUtils/man/isotopologues.html),
    which checks mass differences and intensity ratios, and
2.  **adducts and fragments**, by neutral-mass voting with a rule list:
    a neutral mass is suggested when two or more features of a group can
    be explained as different ions of it.

All supported neutral-mass hypotheses are kept and ranked. The
annotations are meant as **hints** for what to check, so the default
rule list is long and recall is preferred over precision. Results are
written to the
[`featureDefinitions()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
of the `XcmsExperiment`: the `isotopes` and `adduct` columns use
CAMERA’s text syntax, and the grouping itself is the `feature_group`
column written by
[`groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html).
The vignette *How neutral-mass voting works* explains the method step by
step.

### Related packages

- *CAMERA* does grouping, isotope and adduct annotation, but works on
  the legacy `xcmsSet` class. *IonFamilies* re-implements its
  neutral-mass voting for
  [`groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
  groups and keeps all hypotheses, ranked, instead of one per feature.
- *InterpretMSSpectrum* (`findMAIN()`) ranks main-adduct hypotheses for
  a single spectrum. *RAMClustR* uses it after its own clustering.
- *cliqueMS* annotates one sample at a time with its own grouping.
- *MetaboAnnotation* matches features against compound databases. It
  complements the database-free annotation here.

## Installation

``` r

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("stanstrup/IonFamilies")
```

## Grouping features

We use the small example data set that comes with *xcms*: an
`XcmsExperiment` with peak detection, alignment and correspondence
already done. First the features are grouped by retention time, then the
groups are split by the correlation of the features’ abundances across
samples.

``` r

library(xcms)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
library(MsFeatures)
library(IonFamilies)
library(dplyr)

x <- loadXcmsData("xmse")
x <- groupFeatures(x, SimilarRtimeParam(diffRt = 10))
x <- groupFeatures(x, AbundanceSimilarityParam(threshold = 0.7,
                                               transform = log2))
table(lengths(split(featureGroups(x), featureGroups(x))))
#> 
#>   1   2   3   4   6 
#> 184  55  13   3   1
```

With raw data available,
[`EicSimilarityParam()`](https://rdrr.io/pkg/xcms/man/groupFeatures-eic-similarity.html)
would add a third step that compares peak shapes.

## Isotopes

[`annotateIsotopes()`](https://stanstrup.github.io/IonFamilies/reference/annotateIsotopes.md)
searches isotope patterns within each group. The representative
intensity of a feature is its median over the samples given by `subset`
(all samples by default).

``` r

x <- annotateIsotopes(x, IsotopologueParam(ppm = 10, tolerance = 0.005))
fd <- as.data.frame(featureDefinitions(x))
fd %>%
    filter(isotopes != "") %>%
    select(mzmed, rtmed, feature_group, isotopes) %>%
    head()
#>       mzmed    rtmed feature_group   isotopes
#> FT033 300.2 3404.692    FG.002.001    [1][M]+
#> FT034 301.0 2789.217    FG.006.001    [4][M]+
#> FT035 301.2 3404.692    FG.002.001  [1][M+1]+
#> FT037 302.0 2790.972    FG.006.001  [4][M+1]+
#> FT069 329.2 3642.045    FG.089.001   [36][M]+
#> FT070 330.2 3643.682    FG.089.001 [36][M+1]+
```

[`isIsotopologue()`](https://stanstrup.github.io/IonFamilies/reference/isIsotopologue.md)
identifies the heavy isotopologues, for example to remove them before
statistics:

``` r

table(isIsotopologue(fd$isotopes))
#> 
#> FALSE  TRUE 
#>   307    44
```

## Adducts and fragments

[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
uses the rules in
[`defaultAdductRules()`](https://stanstrup.github.io/IonFamilies/reference/adductRules.md),
which come from the *commonMZ* package and are expanded to one to four
molecules per ion.

``` r

rules <- defaultAdductRules("positive")
nrow(rules)
#> [1] 159
head(rules)
#>        name nmol charge massdiff oidscore quasi  ips
#> 1    [M+H]+    1      1 1.007276        1     1 1.00
#> 2   [2M+H]+    2      1 1.007276        1     0 1.00
#> 3   [3M+H]+    3      1 1.007276        1     0 1.00
#> 4   [4M+H]+    4      1 1.007276        1     0 1.00
#> 5  [M+2H]2+    1      2 2.014552        2     0 0.75
#> 6 [2M+2H]2+    2      2 2.014552        2     0 0.75
```

Every feature of a group is tried as every quasi-molecular ion in the
rules (such as `[M+H]+` or `[M+Na]+`). A neutral mass M becomes a
hypothesis when at least one other feature of the same group is
explained by some rule as an ion of that M. Heavy isotopologues found in
the previous step are left out.

``` r

x <- annotateAdducts(x, AdductVoteParam(polarity = "positive"))
hints <- adductHints(x)
```

[`adductHints()`](https://stanstrup.github.io/IonFamilies/reference/adductHints.md)
returns all hypotheses as one table, with one row per hypothesis and
feature. `hypothesis` ranks the hypotheses within a group by `score`
(the summed rule scores), and `n_ions` is the number of features that
support M.

To see how this works, we look at the group with the most competing
hypotheses:

``` r

group <- hints %>%
    slice_max(hypothesis, n = 1, with_ties = FALSE) %>%
    pull(feature_group)

fd <- featureDefinitions(x) %>%
    as.data.frame() %>%
    mutate(feature = row_number())

fd %>%
    filter(feature_group == group) %>%
    select(mzmed, rtmed, isotopes)
#>       mzmed    rtmed   isotopes
#> FT040 305.1 3504.720           
#> FT093 351.1 3503.953   [15][M]+
#> FT095 353.1 3504.855 [15][M+2]+
```

The group has 3 features. One of them is a heavy isotopologue (`[M+1]+`
or `[M+2]+`) and gets no adduct hypotheses. The hypotheses of the others
are:

``` r

group_hints <- hints %>%
    filter(feature_group == group) %>%
    left_join(select(fd, feature, mzmed), by = "feature") %>%
    arrange(hypothesis, mzmed)

group_hints %>%
    select(hypothesis, M, n_ions, score, mzmed, adduct, rank)
#>   hypothesis        M n_ions score mzmed           adduct rank
#> 1          1 350.0927      2  1.75 305.1     [M+H-HCOOH]+    1
#> 2          1 350.0927      2  1.75 351.1           [M+H]+    1
#> 3          2 282.1108      2  1.50 305.1          [M+Na]+    2
#> 4          2 282.1108      2  1.50 351.1    [M+H+CHOONa]+    2
#> 5          3 368.1033      2  1.25 305.1 [M+H-H2O-HCOOH]+    3
#> 6          3 368.1033      2  1.25 351.1       [M+H-H2O]+    3
```

Each hypothesis explains the same features in a different way. The best
one, M = 350.093, reads them as \[M+H-HCOOH\]+ and \[M+H\]+. The others
are alternatives that the long rule list also allows. This is why the
annotations are hints: they tell you what to check, not what the
compound is. Note that this example data set is small and of low mass
accuracy, so the hypotheses here are weak. The section [Working with
plain vectors](#working-with-plain-vectors) shows a clear case.

### The adduct column

`rank` is the rank of a hypothesis for one feature (1 = most plausible).
The `adduct` column of
[`featureDefinitions()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
shows the `nBest` best hypotheses of each feature as text in CAMERA’s
format, best first. The default is `nBest = 3`.
[`formatAdducts()`](https://stanstrup.github.io/IonFamilies/reference/formatAdducts.md)
rewrites the column with another `nBest` without recomputing anything:

``` r

x1 <- formatAdducts(x, nBest = 1)

fd %>%
    mutate(nBest_1 = featureDefinitions(x1)$adduct) %>%
    filter(feature_group == group) %>%
    select(mzmed, nBest_3 = adduct, nBest_1)
#>       mzmed                                                       nBest_3
#> FT040 305.1 [M+H-HCOOH]+ 350.093 [M+Na]+ 282.111 [M+H-H2O-HCOOH]+ 368.103
#> FT093 351.1       [M+H]+ 350.093 [M+H+CHOONa]+ 282.111 [M+H-H2O]+ 368.103
#> FT095 353.1                                                              
#>                    nBest_1
#> FT040 [M+H-HCOOH]+ 350.093
#> FT093       [M+H]+ 350.093
#> FT095
```

### One hypothesis per feature

[`bestHypotheses()`](https://stanstrup.github.io/IonFamilies/reference/bestHypotheses.md)
reduces the hints to at most one neutral mass per feature, as CAMERA’s
final result did. Hypotheses are accepted from the best down, and one is
skipped if any of its features is already explained:

``` r

bestHypotheses(hints) %>%
    filter(feature_group == group)
#>   feature feature_group hypothesis rank        M       adduct n_ions score
#> 1      40    FG.030.001          1    1 350.0927 [M+H-HCOOH]+      2  1.75
#> 2      93    FG.030.001          1    1 350.0927       [M+H]+      2  1.75
#>   intensity
#> 1  822876.7
#> 2  822876.7
```

### Stored parameters

The parameters used are stored in the object:

``` r

S4Vectors::metadata(x)$IonFamilies
#> $isotopes
#> Object of class "IsotopologueParam"
#>  ppm: 10 
#>  tolerance: 0.005 
#>  polarity: positive 
#>  value: into 
#>  subset: all samples 
#> 
#> $adducts
#> Object of class "AdductVoteParam"
#>  rules: 159 rules, 5 quasi-molecular
#>  ppm: 10 
#>  mzabs: 0.005 
#>  polarity: positive 
#>  value: into 
#>  subset: all samples 
#>  nBest: 3
```

## Working with plain vectors

All functions also work without *xcms*. Here are the ions of
2-furoylglycine (M = 169.0375), a marker of coffee intake, from one
feature group:

``` r

mz <- c(95.0124, 124.0392, 152.0342, 170.0445, 171.0496, 192.0275,
        377.0300, 393.0019)
int <- c(71.6, 81.4, 49.8, 100, 4.33, 11.3, 1.76, 1.82)
grp <- rep("FG.1", length(mz))

iso <- annotateIsotopes(mz, IsotopologueParam(), intensity = int,
                        groups = grp)
hints <- annotateAdducts(mz, AdductVoteParam(), intensity = int,
                         groups = grp, exclude = isIsotopologue(iso))
data.frame(mz, isotopes = iso,
           adduct = formatAdducts(hints, nBest = 2, nFeatures = length(mz)))
#>         mz  isotopes
#> 1  95.0124          
#> 2 124.0392          
#> 3 152.0342          
#> 4 170.0445   [1][M]+
#> 5 171.0496 [1][M+1]+
#> 6 192.0275          
#> 7 377.0300          
#> 8 393.0019          
#>                                                                     adduct
#> 1 [M+H - O=C=NCH3 (methyl isocyanate)]+ 151.027 [M+H-NH3-H2O-H2O]+ 147.055
#> 2                                   [M+H-HCOOH]+ 169.037 [M+H-CO]+ 151.027
#> 3                                        [M+H-H2O]+ 169.037 [M+H]+ 151.027
#> 4                                           [M+H]+ 169.037 [M+Na]+ 147.055
#> 5                                                                         
#> 6                                       [M+Na]+ 169.037 [M+2Na-H]+ 147.055
#> 7                                           [2M+K]+ 169.037 [4M+H]+ 94.005
#> 8                                                           [M+K]+ 354.041
```

Here the best hypothesis is the right one. Five ions agree on M =
169.037: the protonated molecule, the sodium adduct, the potassium
adduct of the dimer, and two in-source fragments (losses of water and of
formic acid). The feature at m/z 171.050 is the `[M+1]+` isotope of the
protonated molecule and is therefore left out:

``` r

hints %>%
    filter(hypothesis == 1) %>%
    mutate(mz = mz[feature], M = round(M, 4)) %>%
    select(mz, adduct, M)
#>         mz       adduct        M
#> 1 124.0392 [M+H-HCOOH]+ 169.0372
#> 2 152.0342   [M+H-H2O]+ 169.0372
#> 3 170.0445       [M+H]+ 169.0372
#> 4 192.0275      [M+Na]+ 169.0372
#> 5 377.0300      [2M+K]+ 169.0372
```

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.1         IonFamilies_0.99.2  MsFeatures_1.20.0  
#> [4] xcms_4.10.1         BiocParallel_1.46.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] DBI_1.3.0                   rlang_1.3.0                
#>  [3] magrittr_2.0.5              clue_0.3-68                
#>  [5] MassSpecWavelet_1.78.2      otel_0.3.0                 
#>  [7] matrixStats_1.5.0           compiler_4.6.2             
#>  [9] PTMods_1.0.0                vctrs_0.7.3                
#> [11] reshape2_1.4.5              stringr_1.6.0              
#> [13] ProtGenerics_1.44.0         pkgconfig_2.0.3            
#> [15] MetaboCoreUtils_1.20.1      crayon_1.5.3               
#> [17] fastmap_1.2.0               XVector_0.52.0             
#> [19] rmarkdown_2.32              preprocessCore_1.74.0      
#> [21] purrr_1.2.2                 xfun_0.62                  
#> [23] MultiAssayExperiment_1.38.0 jsonlite_2.0.0             
#> [25] progress_1.2.3              DelayedArray_0.38.2        
#> [27] parallel_4.6.2              prettyunits_1.2.0          
#> [29] cluster_2.1.8.3             R6_2.6.1                   
#> [31] stringi_1.8.9               RColorBrewer_1.1-3         
#> [33] limma_3.68.5                GenomicRanges_1.64.0       
#> [35] Rcpp_1.1.2                  Seqinfo_1.2.0              
#> [37] SummarizedExperiment_1.42.0 iterators_1.0.14           
#> [39] knitr_1.52                  IRanges_2.46.0             
#> [41] Matrix_1.7-6                igraph_2.3.4               
#> [43] tidyselect_1.2.1            abind_1.4-8                
#> [45] yaml_2.3.12                 doParallel_1.0.17          
#> [47] codetools_0.2-20            affy_1.90.0                
#> [49] lattice_0.23-1              tibble_3.3.1               
#> [51] plyr_1.8.9                  Biobase_2.72.0             
#> [53] withr_3.0.3                 S7_0.2.2                   
#> [55] evaluate_1.0.5              Spectra_1.22.2             
#> [57] pillar_1.11.1               affyio_1.82.0              
#> [59] BiocManager_1.30.27         MatrixGenerics_1.24.0      
#> [61] foreach_1.5.2               stats4_4.6.2               
#> [63] MSnbase_2.38.1              MALDIquant_1.22.3          
#> [65] ncdf4_1.24                  generics_0.1.4             
#> [67] S4Vectors_0.50.3            hms_1.1.4                  
#> [69] ggplot2_4.0.3               scales_1.4.0               
#> [71] MsExperiment_1.14.0         glue_1.8.1                 
#> [73] lazyeval_0.2.3              tools_4.6.2                
#> [75] mzID_1.50.0                 data.table_1.18.6.1        
#> [77] QFeatures_1.22.0            vsn_3.80.0                 
#> [79] mzR_2.46.0                  fs_2.1.0                   
#> [81] XML_3.99-0.25               grid_4.6.2                 
#> [83] impute_1.86.0               tidyr_1.3.2                
#> [85] MsCoreUtils_1.24.0          PSMatch_1.16.0             
#> [87] cli_3.6.6                   S4Arrays_1.12.1            
#> [89] AnnotationFilter_1.36.0     pcaMethods_2.4.0           
#> [91] gtable_0.3.6                digest_0.6.39              
#> [93] BiocGenerics_0.58.1         SparseArray_1.12.3         
#> [95] farver_2.1.2                htmltools_0.5.9            
#> [97] lifecycle_1.0.5             statmod_1.5.2              
#> [99] MASS_7.3-66
```
