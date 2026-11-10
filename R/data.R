#' Default adduct and fragment rules
#'
#' `adductRules` is a list with two `data.frame`s, `positive` and `negative`,
#' with the adduct and in-source fragment rules from the 'commonMZ' package
#' (`commonMZ::MZ_CAMERA()`), in CAMERA's rule format and expanded to one to
#' four molecules per ion. The list is long on purpose: annotations are
#' meant as hints, so a missed true relationship is worse than a coincidental
#' one. The ranking of [annotateAdducts()] puts the best supported hypothesis
#' first.
#'
#' `defaultAdductRules()` returns the table for one polarity.
#'
#' @format A list of two `data.frame`s with the columns:
#' \describe{
#'   \item{name}{rule name, e.g. `"[M+H]+"`.}
#'   \item{nmol}{number of molecules in the ion.}
#'   \item{charge}{charge of the ion.}
#'   \item{massdiff}{mass added to `nmol` neutral molecules to form the ion.}
#'   \item{oidscore}{CAMERA's rule identifier score.}
#'   \item{quasi}{1 for quasi-molecular ions, which can anchor a
#'     hypothesis.}
#'   \item{ips}{rule score; higher means more common.}
#' }
#'
#' @source Generated from 'commonMZ' (MIT licence, by the same author) with
#'     `inst/scripts/make-adductRules.R`.
"adductRules"

#' @rdname adductRules
#'
#' @param polarity `character(1)`, `"positive"` or `"negative"`.
#'
#' @return `defaultAdductRules()`: a `data.frame`.
#'
#' @export
#'
#' @examples
#' head(defaultAdductRules("positive"))
defaultAdductRules <- function(polarity = c("positive", "negative")) {
    polarity <- match.arg(polarity)
    env <- new.env(parent = emptyenv())
    utils::data("adductRules", package = "IonFamilies", envir = env)
    env$adductRules[[polarity]]
}
