## Parameter classes ----------------------------------------------------------

#' Parameters for isotope annotation
#'
#' `IsotopologueParam()` creates the parameter object for [annotateIsotopes()].
#' Isotope patterns are searched within each feature group with
#' [MetaboCoreUtils::isotopologues()], which checks both the mass differences
#' and the intensity ratios of the candidate peaks.
#'
#' @param ppm `numeric(1)` relative m/z tolerance (in ppm) passed to
#'     [MetaboCoreUtils::isotopologues()].
#' @param tolerance `numeric(1)` absolute m/z tolerance passed to
#'     [MetaboCoreUtils::isotopologues()].
#' @param polarity `character(1)`, `"positive"` or `"negative"`. Only used for
#'     the `+`/`-` suffix of the labels.
#' @param subset `integer` with the indices of the samples used to compute
#'     the representative intensity of each feature (the median over these
#'     samples). An empty vector (the default) uses all samples. Only used by
#'     the `XcmsExperiment` method.
#' @param value `character(1)` the type of feature value used for the
#'     intensities, see [xcms::featureValues()]. Only used by the
#'     `XcmsExperiment` method.
#'
#' @return An `IsotopologueParam` object.
#'
#' @seealso [annotateIsotopes()], [AdductVoteParam()]
#'
#' @export
#'
#' @examples
#' IsotopologueParam()
#' IsotopologueParam(ppm = 5, polarity = "negative")
IsotopologueParam <- function(ppm = 10, tolerance = 0.005,
                              polarity = c("positive", "negative"),
                              subset = integer(), value = "into") {
    polarity <- match.arg(polarity)
    new("IsotopologueParam", ppm = ppm, tolerance = tolerance,
        polarity = polarity, subset = as.integer(subset), value = value)
}

#' @rdname IsotopologueParam
#' @exportClass IsotopologueParam
setClass("IsotopologueParam",
         slots = c(ppm = "numeric", tolerance = "numeric",
                   polarity = "character", subset = "integer",
                   value = "character"),
         prototype = prototype(ppm = 10, tolerance = 0.005,
                               polarity = "positive", subset = integer(),
                               value = "into"))

setValidity("IsotopologueParam", function(object) {
    msg <- c(.check_tolerance(object@ppm, "ppm"),
             .check_tolerance(object@tolerance, "tolerance"),
             .check_polarity(object@polarity),
             .check_value(object@value))
    if (length(msg)) msg else TRUE
})

#' Parameters for adduct and fragment annotation
#'
#' `AdductVoteParam()` creates the parameter object for [annotateAdducts()].
#'
#' Within each feature group, every feature is tried as every
#' quasi-molecular ion of the rule list (rules with `quasi == 1`, such as
#' `[M+H]+`). The neutral mass this implies is a hypothesis if at least one
#' other feature of the group is explained by some rule as an ion of the
#' same neutral mass. All hypotheses are kept and ranked by the summed rule
#' score (`ips`) and then by the summed intensity of the features they
#' explain.
#'
#' @param rules `data.frame` with the adduct and fragment rules in CAMERA's
#'     format, with columns `name`, `nmol`, `charge`, `massdiff`, `quasi` and
#'     `ips`. `NULL` (the default) uses [defaultAdductRules()] for
#'     `polarity`.
#' @param ppm `numeric(1)` relative tolerance (in ppm) for two neutral masses
#'     to be considered equal.
#' @param mzabs `numeric(1)` absolute tolerance added to `ppm`.
#' @param polarity `character(1)`, `"positive"` or `"negative"`. Selects the
#'     default rules.
#' @param subset `integer` with the indices of the samples used to compute
#'     the representative intensity of each feature. An empty vector (the
#'     default) uses all samples. Only used by the `XcmsExperiment` method.
#' @param value `character(1)` the type of feature value used for the
#'     intensities, see [xcms::featureValues()]. Only used by the
#'     `XcmsExperiment` method.
#' @param nBest `numeric(1)` how many hypotheses per feature are written to
#'     the `adduct` text column. `Inf` writes all of them. All hypotheses are
#'     stored regardless; see [adductHints()] and [formatAdducts()].
#'
#' @return An `AdductVoteParam` object.
#'
#' @seealso [annotateAdducts()], [defaultAdductRules()], [IsotopologueParam()]
#'
#' @export
#'
#' @examples
#' AdductVoteParam()
#' AdductVoteParam(polarity = "negative", nBest = Inf)
AdductVoteParam <- function(rules = NULL, ppm = 10, mzabs = 0.005,
                            polarity = c("positive", "negative"),
                            subset = integer(), value = "into", nBest = 3) {
    polarity <- match.arg(polarity)
    if (is.null(rules))
        rules <- defaultAdductRules(polarity)
    new("AdductVoteParam", rules = as.data.frame(rules), ppm = ppm,
        mzabs = mzabs, polarity = polarity, subset = as.integer(subset),
        value = value, nBest = nBest)
}

#' @rdname AdductVoteParam
#' @exportClass AdductVoteParam
setClass("AdductVoteParam",
         slots = c(rules = "data.frame", ppm = "numeric", mzabs = "numeric",
                   polarity = "character", subset = "integer",
                   value = "character", nBest = "numeric"),
         prototype = prototype(rules = data.frame(name = character(),
                                                  nmol = integer(),
                                                  charge = integer(),
                                                  massdiff = numeric(),
                                                  quasi = integer(),
                                                  ips = numeric()),
                               ppm = 10, mzabs = 0.005,
                               polarity = "positive", subset = integer(),
                               value = "into", nBest = 3))

setValidity("AdductVoteParam", function(object) {
    msg <- c(.check_rules(object@rules),
             .check_tolerance(object@ppm, "ppm"),
             .check_tolerance(object@mzabs, "mzabs"),
             .check_polarity(object@polarity),
             .check_value(object@value),
             .check_n_best(object@nBest))
    if (length(msg)) msg else TRUE
})

## Validity helpers -----------------------------------------------------------

.RULE_COLUMNS <- c("name", "nmol", "charge", "massdiff", "quasi", "ips")

.check_tolerance <- function(x, name) {
    if (length(x) != 1L || is.na(x) || x < 0)
        paste0("'", name, "' has to be a single non-negative number")
    else character()
}

.check_polarity <- function(x) {
    if (length(x) != 1L || !x %in% c("positive", "negative"))
        "'polarity' has to be \"positive\" or \"negative\""
    else character()
}

.check_value <- function(x) {
    if (length(x) != 1L || is.na(x))
        "'value' has to be a single character string"
    else character()
}

.check_n_best <- function(x) {
    if (length(x) != 1L || is.na(x) || x < 1)
        "'nBest' has to be a single number >= 1 (Inf for all)"
    else character()
}

.check_rules <- function(x) {
    missing_cols <- setdiff(.RULE_COLUMNS, colnames(x))
    if (length(missing_cols))
        paste0("'rules' lacks the column(s) ",
               paste0("'", missing_cols, "'", collapse = ", "))
    else character()
}
