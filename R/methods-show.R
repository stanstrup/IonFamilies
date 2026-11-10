## show() methods -------------------------------------------------------------

#' @rdname IsotopologueParam
#' @param object an `IsotopologueParam` object.
setMethod("show", "IsotopologueParam", function(object) {
    cat("Object of class \"IsotopologueParam\"\n")
    cat(" ppm:", object@ppm, "\n")
    cat(" tolerance:", object@tolerance, "\n")
    cat(" polarity:", object@polarity, "\n")
    cat(" value:", object@value, "\n")
    cat(" subset:", .describe_subset(object@subset), "\n")
})

#' @rdname AdductVoteParam
#' @param object an `AdductVoteParam` object.
setMethod("show", "AdductVoteParam", function(object) {
    cat("Object of class \"AdductVoteParam\"\n")
    cat(" rules:", nrow(object@rules), "rules,",
        sum(object@rules$quasi == 1L), "quasi-molecular\n")
    cat(" ppm:", object@ppm, "\n")
    cat(" mzabs:", object@mzabs, "\n")
    cat(" polarity:", object@polarity, "\n")
    cat(" value:", object@value, "\n")
    cat(" subset:", .describe_subset(object@subset), "\n")
    cat(" nBest:", object@nBest, "\n")
})

.describe_subset <- function(x) {
    if (!length(x))
        return("all samples")
    paste(length(x), "samples")
}
