## Methods for plain vectors --------------------------------------------------

#' @rdname annotateIsotopes
setMethod("annotateIsotopes", signature("numeric", "IsotopologueParam"),
          function(object, param, intensity, groups, ...) {
              .check_lengths(object, intensity, groups)
              .annotate_isotopes(object, intensity, groups,
                                 ppm = param@ppm,
                                 tolerance = param@tolerance,
                                 polarity = param@polarity)
          })

#' @rdname annotateAdducts
setMethod("annotateAdducts", signature("numeric", "AdductVoteParam"),
          function(object, param, intensity, groups,
                   exclude = rep(FALSE, length(object)), ...) {
              .check_lengths(object, intensity, groups)
              if (length(exclude) != length(object))
                  stop("'exclude' must have the same length as 'object'")
              .annotate_adducts(object, intensity, groups,
                                rules = param@rules, exclude = exclude,
                                ppm = param@ppm, mzabs = param@mzabs)
          })

#' @rdname formatAdducts
setMethod("formatAdducts", "data.frame",
          function(object, nBest = 3, nFeatures = max(object$feature, 0L),
                   ...) {
              msg <- .check_n_best(nBest)
              if (length(msg))
                  stop(msg)
              .format_hints(object, nFeatures, nBest)
          })

.check_lengths <- function(mz, intensity, groups) {
    if (length(intensity) != length(mz) || length(groups) != length(mz))
        stop("'intensity' and 'groups' must have the same length as ",
             "'object' (the m/z values)")
    invisible(TRUE)
}
