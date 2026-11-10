## Methods for XcmsExperiment -------------------------------------------------

#' @rdname annotateIsotopes
#' @importClassesFrom xcms XcmsExperiment
setMethod("annotateIsotopes", signature("XcmsExperiment", "IsotopologueParam"),
          function(object, param, ...) {
              .check_feature_groups(object)
              fd <- featureDefinitions(object)
              intensity <- .representative_intensity(object, param@value,
                                                     param@subset)
              fd$isotopes <- .annotate_isotopes(fd$mzmed, intensity,
                                                featureGroups(object),
                                                ppm = param@ppm,
                                                tolerance = param@tolerance,
                                                polarity = param@polarity)
              featureDefinitions(object) <- fd
              .store_param(object, "isotopes", param)
          })

#' @rdname annotateAdducts
setMethod("annotateAdducts", signature("XcmsExperiment", "AdductVoteParam"),
          function(object, param, ...) {
              .check_feature_groups(object)
              fd <- featureDefinitions(object)
              groups <- featureGroups(object)
              intensity <- .representative_intensity(object, param@value,
                                                     param@subset)
              if ("isotopes" %in% colnames(fd)) {
                  exclude <- isIsotopologue(fd$isotopes)
              } else {
                  message("No isotope annotation found; consider running ",
                          "'annotateIsotopes()' first so that heavy ",
                          "isotopologues are not used as adducts.")
                  exclude <- rep(FALSE, nrow(fd))
              }
              hints <- .annotate_adducts(fd$mzmed, intensity, groups,
                                         rules = param@rules,
                                         exclude = exclude, ppm = param@ppm,
                                         mzabs = param@mzabs)
              fd$pcgroup <- .pcgroup_numbers(groups, intensity)
              fd$adduct <- .format_hints(hints, nrow(fd), param@nBest)
              fd$adduct_hints <- I(.hints_to_list(hints, nrow(fd)))
              featureDefinitions(object) <- fd
              .store_param(object, "adducts", param)
          })

#' @rdname adductHints
setMethod("adductHints", "XcmsExperiment", function(object, ...) {
    fd <- featureDefinitions(object)
    if (!"adduct_hints" %in% colnames(fd))
        stop("No adduct annotation found. Run 'annotateAdducts()' first.")
    .list_to_hints(fd$adduct_hints)
})

#' @rdname formatAdducts
setMethod("formatAdducts", "XcmsExperiment",
          function(object, nBest = 3, ...) {
              msg <- .check_n_best(nBest)
              if (length(msg))
                  stop(msg)
              hints <- adductHints(object)
              fd <- featureDefinitions(object)
              fd$adduct <- .format_hints(hints, nrow(fd), nBest)
              featureDefinitions(object) <- fd
              param <- metadata(object)$IonFamilies$adducts
              if (is(param, "AdductVoteParam")) {
                  param@nBest <- nBest
                  object <- .store_param(object, "adducts", param)
              }
              object
          })

## Internal helpers -----------------------------------------------------------

.check_feature_groups <- function(x) {
    fd <- featureDefinitions(x)
    if (!nrow(fd))
        stop("No features present. Run 'groupChromPeaks()' first.")
    if (!"feature_group" %in% colnames(fd))
        stop("No feature groups present. Run 'groupFeatures()' first.")
    invisible(TRUE)
}

## Median feature value over the selected samples; NA becomes 0
.representative_intensity <- function(x, value, subset) {
    vals <- featureValues(x, value = value)
    if (length(subset))
        vals <- vals[, subset, drop = FALSE]
    res <- apply(vals, 1L, median, na.rm = TRUE)
    res[is.na(res)] <- 0
    unname(res)
}

.store_param <- function(x, step, param) {
    ion_families <- metadata(x)$IonFamilies
    if (is.null(ion_families))
        ion_families <- list()
    ion_families[[step]] <- param
    metadata(x)$IonFamilies <- ion_families
    x
}
