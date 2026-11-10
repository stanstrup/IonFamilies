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
of the `XcmsExperiment` in the format CAMERA used (`pcgroup`,
`isotopes`, `adduct`), so existing downstream code keeps working. The
vignette *How neutral-mass voting works* explains the method step by
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
library(MsFeatures)
library(IonFamilies)

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
fd <- featureDefinitions(x)
head(fd[fd$isotopes != "", c("mzmed", "rtmed", "feature_group", "isotopes")])
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

`nBest` sets how many hypotheses per feature are written to the `adduct`
column. All of them are stored in any case.

``` r

x <- annotateAdducts(x, AdductVoteParam(polarity = "positive", nBest = 2))
fd <- featureDefinitions(x)
head(fd[fd$adduct != "", c("mzmed", "pcgroup", "isotopes", "adduct")])
#>       mzmed pcgroup isotopes                                 adduct
#> FT002 205.0      55                            [M+H-C3H9N]+ 263.066
#> FT003 206.0     146                            [M+H-C4H6O]+ 275.037
#> FT004 207.1      80                      [M+H-S-NH3-HCOOH]+ 301.103
#> FT021 281.2      41          [M+H-HCOOH]+ 326.203 [M+H-CO]+ 308.193
#> FT022 281.1      55                                [M+NH4]+ 263.066
#> FT027 284.1      80                              [M+H-H2O]+ 301.103
```

The full table of hypotheses, one row per hypothesis and feature, is
returned by
[`adductHints()`](https://stanstrup.github.io/IonFamilies/reference/adductHints.md).
`rank` is the rank of the hypothesis for that feature and `n_ions` the
number of features that support the neutral mass.

``` r

hints <- adductHints(x)
head(hints)
#>   feature      group hypothesis rank        M             adduct n_ions score
#> 1       2 FG.005.002          1    1 263.0662       [M+H-C3H9N]+      2  1.50
#> 2       3 FG.005.001          1    1 275.0368       [M+H-C4H6O]+      2  1.50
#> 3       4 FG.078.001          1    1 301.1033 [M+H-S-NH3-HCOOH]+      2  1.25
#> 4      21 FG.001.001          1    1 326.2033       [M+H-HCOOH]+      3  1.75
#> 5      21 FG.001.001          2    2 308.1927          [M+H-CO]+      2  1.75
#> 6      22 FG.005.002          1    1 263.0662           [M+NH4]+      2  1.50
#>   intensity
#> 1 1948834.6
#> 2  501911.2
#> 3 1093366.1
#> 4 2844075.4
#> 5 2341545.7
#> 6 1948834.6
```

The text column can be rewritten with another `nBest` without
recomputing:

``` r

x <- formatAdducts(x, nBest = 1)
head(featureDefinitions(x)$adduct[featureDefinitions(x)$adduct != ""])
#> [1] "[M+H-C3H9N]+ 263.066"       "[M+H-C4H6O]+ 275.037"      
#> [3] "[M+H-S-NH3-HCOOH]+ 301.103" "[M+H-HCOOH]+ 326.203"      
#> [5] "[M+NH4]+ 263.066"           "[M+H-H2O]+ 301.103"
```

[`bestHypotheses()`](https://stanstrup.github.io/IonFamilies/reference/bestHypotheses.md)
reduces the hints to at most one neutral mass per feature, as CAMERA’s
final result did:

``` r

head(bestHypotheses(hints))
#>   feature      group hypothesis rank        M       adduct n_ions score
#> 1      21 FG.001.001          1    1 326.2033 [M+H-HCOOH]+      3  1.75
#> 2      46 FG.001.001          1    1 326.2033   [M+H-H2O]+      3  1.75
#> 3      76 FG.001.001          1    1 326.2033  [2M+H+Na]2+      3  1.75
#> 4       3 FG.005.001          1    1 275.0368 [M+H-C4H6O]+      2  1.50
#> 5      54 FG.005.001          1    1 275.0368       [M+K]+      2  1.50
#> 6       2 FG.005.002          1    1 263.0662 [M+H-C3H9N]+      2  1.50
#>   intensity
#> 1 2844075.4
#> 2 2844075.4
#> 3 2844075.4
#> 4  501911.2
#> 5  501911.2
#> 6 1948834.6
```

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
#>  nBest: 1
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
#> [1] IonFamilies_0.99.0  MsFeatures_1.20.0   xcms_4.10.1        
#> [4] BiocParallel_1.46.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1            dplyr_1.2.1                
#>  [3] farver_2.1.2                S7_0.2.2                   
#>  [5] fastmap_1.2.0               lazyeval_0.2.3             
#>  [7] MassSpecWavelet_1.78.2      XML_3.99-0.25              
#>  [9] digest_0.6.39               lifecycle_1.0.5            
#> [11] cluster_2.1.8.3             ProtGenerics_1.44.0        
#> [13] PTMods_1.0.0                statmod_1.5.2              
#> [15] magrittr_2.0.5              compiler_4.6.2             
#> [17] progress_1.2.3              rlang_1.3.0                
#> [19] tools_4.6.2                 igraph_2.3.4               
#> [21] yaml_2.3.12                 data.table_1.18.6.1        
#> [23] knitr_1.52                  prettyunits_1.2.0          
#> [25] S4Arrays_1.12.1             DelayedArray_0.38.2        
#> [27] plyr_1.8.9                  RColorBrewer_1.1-3         
#> [29] abind_1.4-8                 purrr_1.2.2                
#> [31] BiocGenerics_0.58.1         grid_4.6.2                 
#> [33] stats4_4.6.2                preprocessCore_1.74.0      
#> [35] ggplot2_4.0.3               iterators_1.0.14           
#> [37] scales_1.4.0                MASS_7.3-66                
#> [39] MultiAssayExperiment_1.38.0 SummarizedExperiment_1.42.0
#> [41] cli_3.6.6                   crayon_1.5.3               
#> [43] mzR_2.46.0                  rmarkdown_2.32             
#> [45] generics_0.1.4              otel_0.3.0                 
#> [47] reshape2_1.4.5              ncdf4_1.24                 
#> [49] DBI_1.3.0                   affy_1.90.0                
#> [51] stringr_1.6.0               parallel_4.6.2             
#> [53] impute_1.86.0               BiocManager_1.30.27        
#> [55] AnnotationFilter_1.36.0     XVector_0.52.0             
#> [57] vsn_3.80.0                  matrixStats_1.5.0          
#> [59] vctrs_0.7.3                 Matrix_1.7-6               
#> [61] jsonlite_2.0.0              hms_1.1.4                  
#> [63] IRanges_2.46.0              S4Vectors_0.50.3           
#> [65] MsExperiment_1.14.0         MALDIquant_1.22.3          
#> [67] clue_0.3-68                 foreach_1.5.2              
#> [69] limma_3.68.5                tidyr_1.3.2                
#> [71] affyio_1.82.0               glue_1.8.1                 
#> [73] MSnbase_2.38.1              codetools_0.2-20           
#> [75] QFeatures_1.22.0            Spectra_1.22.2             
#> [77] stringi_1.8.9               gtable_0.3.6               
#> [79] GenomicRanges_1.64.0        mzID_1.50.0                
#> [81] tibble_3.3.1                pillar_1.11.1              
#> [83] pcaMethods_2.4.0            htmltools_0.5.9            
#> [85] Seqinfo_1.2.0               R6_2.6.1                   
#> [87] doParallel_1.0.17           evaluate_1.0.5             
#> [89] lattice_0.23-1              Biobase_2.72.0             
#> [91] MetaboCoreUtils_1.20.1      Rcpp_1.1.2                 
#> [93] PSMatch_1.16.0              SparseArray_1.12.3         
#> [95] xfun_0.62                   MsCoreUtils_1.24.0         
#> [97] MatrixGenerics_1.24.0       fs_2.1.0                   
#> [99] pkgconfig_2.0.3
```
