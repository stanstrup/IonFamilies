## Isotope annotation ---------------------------------------------------------

## Mass difference between 13C and 12C
.C13_STEP <- 1.003355

#' Is an isotope label a heavy isotopologue?
#'
#' `isIsotopologue()` returns `TRUE` for labels of heavy isotopologues
#' (`"[3][M+1]+"`, `"[3][M+2]-"`) and `FALSE` for monoisotopic peaks
#' (`"[3][M]+"`), adduct labels such as `"[M+H]+"`, and `""`.
#'
#' @param x `character` with isotope labels as written by
#'     [annotateIsotopes()] or CAMERA.
#'
#' @return `logical` of the same length as `x`.
#'
#' @export
#'
#' @examples
#' isIsotopologue(c("[1][M]+", "[1][M+1]+", "[M+H]+", ""))
isIsotopologue <- function(x) {
    grepl("\\[M\\+[0-9]+\\]", x)
}

## Isotope patterns in one group. `idx` are the indices of the group's
## features in the full vectors; returns a list of index vectors, each one
## pattern with the monoisotopic peak first.
.isotope_patterns_group <- function(idx, mz, intensity, ppm, tolerance) {
    if (length(idx) < 2L)
        return(list())
    ord <- idx[order(mz[idx])]
    peaks <- cbind(mz = mz[ord], intensity = intensity[ord])
    patterns <- isotopologues(peaks, ppm = ppm, tolerance = tolerance)
    lapply(patterns, .subset_index, ord)
}

.subset_index <- function(i, idx) {
    idx[i]
}

## CAMERA-style labels for one pattern with running number `id`
.pattern_labels <- function(members, id, mz, suffix) {
    step <- round((mz[members] - mz[members[1L]]) / .C13_STEP)
    shift <- ifelse(step == 0, "", paste0("+", step))
    paste0("[", id, "][M", shift, "]", suffix)
}

## Workhorse of the numeric method
.annotate_isotopes <- function(mz, intensity, groups, ppm, tolerance,
                               polarity) {
    labels <- character(length(mz))
    by_group <- split(seq_along(mz), groups)
    patterns <- lapply(by_group, .isotope_patterns_group, mz = mz,
                       intensity = intensity, ppm = ppm,
                       tolerance = tolerance)
    patterns <- unlist(unname(patterns), recursive = FALSE)
    if (!length(patterns))
        return(labels)
    suffix <- if (polarity == "positive") "+" else "-"
    pattern_labels <- mapply(.pattern_labels, patterns,
                             seq_along(patterns),
                             MoreArgs = list(mz = mz, suffix = suffix),
                             SIMPLIFY = FALSE)
    labels[unlist(patterns)] <- unlist(pattern_labels)
    labels
}
