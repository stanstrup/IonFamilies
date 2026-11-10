## The ions of 2-furoylglycine (M = 169.0375) from one feature group of the
## coffee intervention study used in the "Introduction to Nutritional
## Metabolomics" course (positive mode). Intensities are relative medians.
furo_mz <- c(95.0124, 124.0392, 152.0342, 170.0445, 171.0496, 192.0275,
             377.0300, 393.0019)
furo_int <- c(71.6, 81.4, 49.8, 100, 4.33, 11.3, 1.76, 1.82)
furo_grp <- rep("FG.1", length(furo_mz))
furo_M <- 169.0375

## Small grouped XcmsExperiment from xcms' example data. loadXcmsData() uses
## data(), which needs xcms attached.
suppressPackageStartupMessages(library(xcms))

grouped_xmse <- function() {
    skip_if_not_installed("faahKO")
    x <- xcms::loadXcmsData("xmse")
    x <- MsFeatures::groupFeatures(x, MsFeatures::SimilarRtimeParam(10))
    MsFeatures::groupFeatures(
        x, MsFeatures::AbundanceSimilarityParam(threshold = 0.7,
                                                transform = log2))
}

n_unique <- function(x) {
    length(unique(x))
}
