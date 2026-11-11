# How neutral-mass voting works

## The example

This vignette walks through the annotation of one feature group step by
step, using the low-level functions that
[`annotateAdducts()`](https://stanstrup.github.io/IonFamilies/reference/annotateAdducts.md)
is built from. The group holds the ions of **2-furoylglycine** (C₇H₇NO₄,
M = 169.0375 Da), a marker of coffee intake, measured in positive mode
in a coffee intervention study. All eight features elute within 0.2 s
and their intensities correlate across the study samples, so
[`groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
put them in one group.

``` r

library(IonFamilies)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
library(dplyr)

mz <- c(95.0124, 124.0392, 152.0342, 170.0445, 171.0496, 192.0275,
        377.0300, 393.0019)
int <- c(71.6, 81.4, 49.8, 100, 4.33, 11.3, 1.76, 1.82)
grp <- rep("FG.1", length(mz))
```

## Step 1: isotopes

Heavy isotopologues are not adducts, so they are found first and left
out of the voting.
[`MetaboCoreUtils::isotopologues()`](https://rdrr.io/pkg/MetaboCoreUtils/man/isotopologues.html)
checks that the mass difference fits an isotope step *and* that the
intensity ratio is plausible for a molecule of that mass.

``` r

iso <- annotateIsotopes(mz, IsotopologueParam(), intensity = int,
                        groups = grp)
data.frame(mz, int, iso)
#>         mz    int       iso
#> 1  95.0124  71.60          
#> 2 124.0392  81.40          
#> 3 152.0342  49.80          
#> 4 170.0445 100.00   [1][M]+
#> 5 171.0496   4.33 [1][M+1]+
#> 6 192.0275  11.30          
#> 7 377.0300   1.76          
#> 8 393.0019   1.82
```

m/z 171.050 is the ¹³C isotopologue of m/z 170.045 (4 % of its
intensity, about what seven carbons give).

## Step 2: implied neutral masses

A rule describes how an ion is formed from `nmol` neutral molecules: its
charge and the mass that is added or lost (`massdiff`). Inverting it
gives the neutral mass a feature would have *if* it were that ion:

``` math
M = \frac{|z| \cdot m/z - \Delta m}{n_{mol}}
```

``` r

rules <- defaultAdductRules("positive")
some <- rules %>%
    filter(name %in% c("[M+H]+", "[M+Na]+", "[M+H-H2O]+", "[M+H-HCOOH]+",
                       "[2M+K]+"))
some %>%
    select(name, nmol, charge, massdiff, quasi, ips)
#>           name nmol charge   massdiff quasi  ips
#> 1       [M+H]+    1      1   1.007276     1 1.00
#> 2      [M+Na]+    1      1  22.989218     1 1.00
#> 3      [2M+K]+    2      1  38.963158     0 0.75
#> 4 [M+H-HCOOH]+    1      1 -44.998203     0 0.75
#> 5   [M+H-H2O]+    1      1 -17.003288     1 0.75
round(impliedMasses(mz, some), 4)
#>        [M+H]+  [M+Na]+  [2M+K]+ [M+H-HCOOH]+ [M+H-H2O]+
#> [1,]  94.0051  72.0232  28.0246     140.0106   112.0157
#> [2,] 123.0319 101.0500  42.5380     169.0374   141.0425
#> [3,] 151.0269 129.0450  56.5355     197.0324   169.0375
#> [4,] 169.0372 147.0553  65.5407     215.0427   187.0478
#> [5,] 170.0423 148.0604  66.0432     216.0478   188.0529
#> [6,] 191.0202 169.0383  76.5322     237.0257   209.0308
#> [7,] 376.0227 354.0408 169.0334     422.0282   394.0333
#> [8,] 391.9946 370.0127 177.0194     438.0001   410.0052
```

Read the table column by column. The value 169.037 appears for m/z
170.045 as `[M+H]+`, for 192.028 as `[M+Na]+`, for 152.034 as
`[M+H-H2O]+`, for 124.039 as `[M+H-HCOOH]+` and for 377.030 as
`[2M+K]+`. Five features agree on one neutral mass.

## Step 3: voting

[`neutralMassVote()`](https://stanstrup.github.io/IonFamilies/reference/neutralMassVote.md)
automates that reading for all rules:

1.  Every feature is tried as every **quasi-molecular** ion (rules with
    `quasi == 1`, such as `[M+H]+` and `[M+Na]+`). Each try is an
    *anchor* and gives a candidate neutral mass.
2.  A candidate becomes a **hypothesis** when at least one *other*
    feature is explained by some rule as an ion of the same mass (within
    `ppm` + `mzabs`).
3.  Hypotheses found from different anchors but with the same members
    are merged.
4.  Hypotheses are **ranked** by the summed rule score (`ips`, higher
    for common ions) and then by the intensity they explain.

``` r

keep <- !isIsotopologue(iso)
votes <- neutralMassVote(mz[keep], int[keep], rules)
votes %>%
    mutate(mz = mz[keep][feature]) %>%
    select(hypothesis, M, n_ions, score, mz, adduct) %>%
    head(12)
#>    hypothesis         M n_ions score       mz
#> 1           1 169.03722      5  4.25 170.0445
#> 2           1 169.03722      5  4.25 192.0275
#> 3           1 169.03722      5  4.25 377.0300
#> 4           1 169.03722      5  4.25 124.0392
#> 5           1 169.03722      5  4.25 152.0342
#> 6           2 151.02692      3  2.25 152.0342
#> 7           2 151.02692      3  2.25 124.0392
#> 8           2 151.02692      3  2.25  95.0124
#> 9           3 147.05528      3  2.25 170.0445
#> 10          3 147.05528      3  2.25 192.0275
#> 11          3 147.05528      3  2.25  95.0124
#> 12          4  94.00512      2  2.00  95.0124
#>                                   adduct
#> 1                                 [M+H]+
#> 2                                [M+Na]+
#> 3                                [2M+K]+
#> 4                           [M+H-HCOOH]+
#> 5                             [M+H-H2O]+
#> 6                                 [M+H]+
#> 7                              [M+H-CO]+
#> 8  [M+H - O=C=NCH3 (methyl isocyanate)]+
#> 9                                [M+Na]+
#> 10                            [M+2Na-H]+
#> 11                    [M+H-NH3-H2O-H2O]+
#> 12                                [M+H]+
```

Hypothesis 1 is furoylglycine with five supporting ions. The later
hypotheses are coincidences: with 159 rules, most pairs of m/z values
can be connected by *some* pair of rules.

## Step 4: keep all hypotheses or only the best?

CAMERA reported, for each feature, only hypotheses that did not conflict
with better ones.
[`bestHypotheses()`](https://stanstrup.github.io/IonFamilies/reference/bestHypotheses.md)
reproduces that choice:

``` r

hints <- annotateAdducts(mz, AdductVoteParam(), intensity = int,
                         groups = grp, exclude = !keep)
bestHypotheses(hints) %>%
    mutate(mz = mz[feature]) %>%
    select(mz, M, adduct, n_ions)
#>         mz        M       adduct n_ions
#> 1 124.0392 169.0372 [M+H-HCOOH]+      5
#> 2 152.0342 169.0372   [M+H-H2O]+      5
#> 3 170.0445 169.0372       [M+H]+      5
#> 4 192.0275 169.0372      [M+Na]+      5
#> 5 377.0300 169.0372      [2M+K]+      5
```

*IonFamilies* keeps **all** hypotheses by default, ranked per feature. A
true annotation that ranks second for a feature is then still visible.
Missing a true relationship costs more than showing an extra suggestion,
because annotations are hints to be checked. `nBest` limits how many are
written to text:

``` r

data.frame(mz, iso,
           best_only = formatAdducts(hints, nBest = 1, nFeatures = length(mz)),
           top_three = formatAdducts(hints, nBest = 3, nFeatures = length(mz)))
#>         mz       iso                                     best_only
#> 1  95.0124           [M+H - O=C=NCH3 (methyl isocyanate)]+ 151.027
#> 2 124.0392                                    [M+H-HCOOH]+ 169.037
#> 3 152.0342                                      [M+H-H2O]+ 169.037
#> 4 170.0445   [1][M]+                                [M+H]+ 169.037
#> 5 171.0496 [1][M+1]+                                              
#> 6 192.0275                                         [M+Na]+ 169.037
#> 7 377.0300                                         [2M+K]+ 169.037
#> 8 393.0019                                          [M+K]+ 354.041
#>                                                                                top_three
#> 1 [M+H - O=C=NCH3 (methyl isocyanate)]+ 151.027 [M+H-NH3-H2O-H2O]+ 147.055 [M+H]+ 94.005
#> 2                        [M+H-HCOOH]+ 169.037 [M+H-CO]+ 151.027 [M+H-H2O-HCOOH]+ 187.048
#> 3                                [M+H-H2O]+ 169.037 [M+H]+ 151.027 [M+H-(H2O)2]+ 187.048
#> 4                                      [M+H]+ 169.037 [M+Na]+ 147.055 [M+H-H2O]+ 187.048
#> 5                                                                                       
#> 6                                 [M+Na]+ 169.037 [M+2Na-H]+ 147.055 [M+Na+K-H]+ 131.081
#> 7                                         [2M+K]+ 169.037 [4M+H]+ 94.005 [M+Na]+ 354.041
#> 8                                                                         [M+K]+ 354.041
```

## Why the grouping matters

Because a long rule list connects almost anything, voting is only
meaningful inside groups of features that are already known to belong
together: same retention time, correlated intensities across samples,
similar peak shapes. That evidence comes from
[`groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html).
*IonFamilies* adds the chemical interpretation on top. Good groups give
good hints.

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
#> [1] dplyr_1.2.1        IonFamilies_0.99.2
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.3.0                   rlang_1.3.0                
#>   [3] magrittr_2.0.5              clue_0.3-68                
#>   [5] MassSpecWavelet_1.78.2      otel_0.3.0                 
#>   [7] matrixStats_1.5.0           compiler_4.6.2             
#>   [9] PTMods_1.0.0                vctrs_0.7.3                
#>  [11] reshape2_1.4.5              stringr_1.6.0              
#>  [13] ProtGenerics_1.44.0         pkgconfig_2.0.3            
#>  [15] MetaboCoreUtils_1.20.1      crayon_1.5.3               
#>  [17] fastmap_1.2.0               XVector_0.52.0             
#>  [19] rmarkdown_2.32              preprocessCore_1.74.0      
#>  [21] xcms_4.10.1                 purrr_1.2.2                
#>  [23] xfun_0.62                   MultiAssayExperiment_1.38.0
#>  [25] jsonlite_2.0.0              progress_1.2.3             
#>  [27] DelayedArray_0.38.2         BiocParallel_1.46.0        
#>  [29] parallel_4.6.2              prettyunits_1.2.0          
#>  [31] cluster_2.1.8.3             R6_2.6.1                   
#>  [33] stringi_1.8.9               RColorBrewer_1.1-3         
#>  [35] limma_3.68.5                GenomicRanges_1.64.0       
#>  [37] Rcpp_1.1.2                  Seqinfo_1.2.0              
#>  [39] SummarizedExperiment_1.42.0 iterators_1.0.14           
#>  [41] knitr_1.52                  IRanges_2.46.0             
#>  [43] Matrix_1.7-6                igraph_2.3.4               
#>  [45] tidyselect_1.2.1            abind_1.4-8                
#>  [47] yaml_2.3.12                 doParallel_1.0.17          
#>  [49] codetools_0.2-20            affy_1.90.0                
#>  [51] lattice_0.23-1              tibble_3.3.1               
#>  [53] plyr_1.8.9                  Biobase_2.72.0             
#>  [55] withr_3.0.3                 S7_0.2.2                   
#>  [57] evaluate_1.0.5              Spectra_1.22.2             
#>  [59] pillar_1.11.1               affyio_1.82.0              
#>  [61] BiocManager_1.30.27         MatrixGenerics_1.24.0      
#>  [63] foreach_1.5.2               stats4_4.6.2               
#>  [65] MSnbase_2.38.1              MALDIquant_1.22.3          
#>  [67] ncdf4_1.24                  generics_0.1.4             
#>  [69] S4Vectors_0.50.3            hms_1.1.4                  
#>  [71] ggplot2_4.0.3               scales_1.4.0               
#>  [73] MsExperiment_1.14.0         glue_1.8.1                 
#>  [75] MsFeatures_1.20.0           lazyeval_0.2.3             
#>  [77] tools_4.6.2                 mzID_1.50.0                
#>  [79] data.table_1.18.6.1         QFeatures_1.22.0           
#>  [81] vsn_3.80.0                  mzR_2.46.0                 
#>  [83] fs_2.1.0                    XML_3.99-0.25              
#>  [85] grid_4.6.2                  impute_1.86.0              
#>  [87] tidyr_1.3.2                 MsCoreUtils_1.24.0         
#>  [89] PSMatch_1.16.0              cli_3.6.6                  
#>  [91] S4Arrays_1.12.1             AnnotationFilter_1.36.0    
#>  [93] pcaMethods_2.4.0            gtable_0.3.6               
#>  [95] digest_0.6.39               BiocGenerics_0.58.1        
#>  [97] SparseArray_1.12.3          farver_2.1.2               
#>  [99] htmltools_0.5.9             lifecycle_1.0.5            
#> [101] statmod_1.5.2               MASS_7.3-66
```
