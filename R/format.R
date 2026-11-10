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

## One data.frame of hints per feature, for the list column
.hints_to_list <- function(hints, nFeatures) {
    cols <- setdiff(.HINT_COLUMNS, "feature")
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

## Back from the list column to one long data.frame
.list_to_hints <- function(x) {
    if (!length(x))
        return(.empty_hints())
    res <- do.call(rbind, mapply(.add_feature_index, x, seq_along(x),
                                 SIMPLIFY = FALSE))
    if (is.null(res) || !nrow(res))
        return(.empty_hints())
    res <- res[, .HINT_COLUMNS, drop = FALSE]
    rownames(res) <- NULL
    res
}

## Integer group numbers; 1 = group containing the most intense feature
.pcgroup_numbers <- function(groups, intensity) {
    levels <- unique(groups[order(-intensity)])
    as.integer(factor(groups, levels = levels))
}
