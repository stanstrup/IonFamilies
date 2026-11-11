## Text output and storage helpers --------------------------------------------

## CAMERA-style text per feature, the nBest best (distinct) hints first
.format_hints <- function(hints, nFeatures, nBest) {
    out <- character(nFeatures)
    if (!nrow(hints))
        return(out)
    hints <- hints[order(hints$feature, hints$rank), , drop = FALSE]
    txt <- sprintf("%s %.3f", hints$adduct, hints$M)
    distinct <- !duplicated(paste(hints$feature, txt))
    hints <- hints[distinct, , drop = FALSE]
    txt <- txt[distinct]
    position <- ave(seq_along(txt), hints$feature, FUN = seq_along)
    shown <- position <= nBest
    joined <- tapply(txt[shown], hints$feature[shown], paste, collapse = " ")
    out[as.integer(names(joined))] <- joined
    out
}

## One data.frame of hints per feature, for the list column. The feature
## index and its feature group follow from the row of featureDefinitions()
## and are not stored.
.hints_to_list <- function(hints, nFeatures) {
    cols <- setdiff(.HINT_COLUMNS, c("feature", "feature_group"))
    per_feature <- split(hints[, cols, drop = FALSE],
                         factor(hints$feature, levels = seq_len(nFeatures)))
    lapply(unname(per_feature), .drop_rownames)
}

.drop_rownames <- function(x) {
    rownames(x) <- NULL
    x
}

.add_feature_index <- function(x, i) {
    x$feature <- rep(i, nrow(x))
    x
}

## Back from the list column to one long data.frame; `groups` are the
## feature groups of all features (featureGroups())
.list_to_hints <- function(x, groups) {
    if (!length(x))
        return(.empty_hints())
    res <- do.call(rbind, mapply(.add_feature_index, x, seq_along(x),
                                 SIMPLIFY = FALSE))
    if (is.null(res) || !nrow(res))
        return(.empty_hints())
    res$feature_group <- groups[res$feature]
    res <- res[, .HINT_COLUMNS, drop = FALSE]
    rownames(res) <- NULL
    res
}
