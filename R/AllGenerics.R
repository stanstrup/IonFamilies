#' Annotate isotopes within feature groups
#'
#' `annotateIsotopes()` labels isotope patterns within each feature group in
#' CAMERA syntax: `"[3][M]+"` is the monoisotopic peak of pattern 3,
#' `"[3][M+1]+"` its first isotopologue, and so on. Features that are not
#' part of a pattern get `""`.
#'
#' @param object For the `numeric` method: the m/z values of the features.
#'     For the `XcmsExperiment` method: an [xcms::XcmsExperiment] with
#'     feature groups (see [MsFeatures::groupFeatures()]).
#' @param param An [IsotopologueParam] object.
#' @param intensity `numeric`, for the `numeric` method only: one
#'     representative intensity per feature, e.g. the median over the study
#'     samples.
#' @param groups for the `numeric` method only: one group label per feature.
#' @param ... ignored.
#'
#' @return For the `numeric` method a `character` vector with one label per
#'     feature. For the `XcmsExperiment` method the object, with the labels in
#'     the `isotopes` column of `featureDefinitions()` and the parameters in
#'     `metadata(object)$IonFamilies$isotopes`.
#'
#' @seealso [IsotopologueParam()], [isIsotopologue()], [annotateAdducts()]
#'
#' @export
#'
#' @examples
#' ## Ions of 2-furoylglycine (M = 169.0375) in one feature group
#' mz <- c(95.0124, 124.0392, 152.0342, 170.0445, 171.0496, 192.0275,
#'         377.0300, 393.0019)
#' int <- c(71.6, 81.4, 49.8, 100, 4.33, 11.3, 1.76, 1.82)
#' annotateIsotopes(mz, IsotopologueParam(), intensity = int,
#'                  groups = rep("FG.1", 8))
setGeneric("annotateIsotopes", function(object, param, ...)
    standardGeneric("annotateIsotopes"))

#' Annotate adducts and fragments within feature groups
#'
#' `annotateAdducts()` finds adducts and in-source fragments by neutral-mass
#' voting within each feature group (see [AdductVoteParam()] for the method).
#' Every supported hypothesis is kept, so a feature can have several ranked
#' alternatives.
#'
#' @param object For the `numeric` method: the m/z values of the features.
#'     For the `XcmsExperiment` method: an [xcms::XcmsExperiment] with
#'     feature groups.
#' @param param An [AdductVoteParam] object.
#' @param intensity `numeric`, for the `numeric` method only: one
#'     representative intensity per feature.
#' @param groups for the `numeric` method only: one group label per feature.
#' @param exclude `logical`, for the `numeric` method only: features to leave
#'     out, normally heavy isotopologues (`isIsotopologue(labels)`). The
#'     `XcmsExperiment` method uses the `isotopes` column if present.
#' @param ... ignored.
#'
#' @return For the `numeric` method a `data.frame` with one row per
#'     hypothesis and feature and the columns `feature` (index),
#'     `feature_group`,
#'     `hypothesis` (rank within the group), `rank` (rank for this feature,
#'     1 = most plausible), `M` (neutral mass), `adduct` (rule name),
#'     `n_ions` (number of features supporting `M`), `score` and
#'     `intensity`. For the `XcmsExperiment` method the object, with the
#'     columns `adduct` (the best `nBest` hints as text) and `adduct_hints`
#'     (a list column with all hypotheses per feature) added to
#'     `featureDefinitions()`, and the parameters in
#'     `metadata(object)$IonFamilies$adducts`.
#'
#' @seealso [AdductVoteParam()], [adductHints()], [formatAdducts()],
#'     [bestHypotheses()], [annotateIsotopes()]
#'
#' @export
#'
#' @examples
#' mz <- c(95.0124, 124.0392, 152.0342, 170.0445, 171.0496, 192.0275,
#'         377.0300, 393.0019)
#' int <- c(71.6, 81.4, 49.8, 100, 4.33, 11.3, 1.76, 1.82)
#' grp <- rep("FG.1", 8)
#' iso <- annotateIsotopes(mz, IsotopologueParam(), intensity = int,
#'                         groups = grp)
#' hints <- annotateAdducts(mz, AdductVoteParam(), intensity = int,
#'                          groups = grp, exclude = isIsotopologue(iso))
#' head(hints)
setGeneric("annotateAdducts", function(object, param, ...)
    standardGeneric("annotateAdducts"))

#' All adduct hypotheses of an annotated object
#'
#' `adductHints()` returns all neutral-mass hypotheses stored by
#' [annotateAdducts()] as one long `data.frame` (one row per hypothesis and
#' feature, columns as described in [annotateAdducts()]).
#'
#' @param object An [xcms::XcmsExperiment] annotated with [annotateAdducts()].
#' @param ... ignored.
#'
#' @return A `data.frame`.
#'
#' @export
#'
#' @examples
#' library(xcms)
#' library(MsFeatures)
#' if (requireNamespace("faahKO", quietly = TRUE)) {
#'     x <- loadXcmsData("xmse")
#'     x <- groupFeatures(x, SimilarRtimeParam(diffRt = 10))
#'     x <- annotateAdducts(x, AdductVoteParam())
#'     head(adductHints(x))
#' }
setGeneric("adductHints", function(object, ...)
    standardGeneric("adductHints"))

#' Write adduct hypotheses as CAMERA-style text
#'
#' `formatAdducts()` writes the `nBest` most plausible hypotheses of each
#' feature as text, best first, in CAMERA's format
#' (`"[M+H]+ 169.037 [M+Na]+ 146.048"`).
#'
#' @param object A `data.frame` of hypotheses as returned by the `numeric`
#'     method of [annotateAdducts()], or an annotated
#'     [xcms::XcmsExperiment].
#' @param nBest `numeric(1)` number of hypotheses per feature; `Inf` for all.
#' @param nFeatures `integer(1)`, for the `data.frame` method only: the total
#'     number of features (length of the result).
#' @param ... ignored.
#'
#' @return For the `data.frame` method a `character` vector with one entry
#'     per feature (`""` without hypotheses). For the `XcmsExperiment` method
#'     the object with a rewritten `adduct` column; nothing is recomputed.
#'
#' @export
#'
#' @examples
#' mz <- c(124.0392, 152.0342, 170.0445, 192.0275, 377.0300)
#' hints <- annotateAdducts(mz, AdductVoteParam(),
#'                          intensity = c(81, 50, 100, 11, 2),
#'                          groups = rep("FG.1", 5))
#' formatAdducts(hints, nBest = 1, nFeatures = 5)
setGeneric("formatAdducts", function(object, ...)
    standardGeneric("formatAdducts"))
