## Neutral-mass voting --------------------------------------------------------

#' Neutral masses implied by m/z values and rules
#'
#' For an ion with m/z `mz` formed by a rule with `nmol` molecules, charge
#' `charge` and mass shift `massdiff`, the neutral mass is
#' `(abs(charge) * mz - massdiff) / nmol`. `impliedMasses()` computes this
#' for every combination of feature and rule. It is the basic operation of
#' rule-based adduct annotation and is exported for teaching and testing.
#'
#' @param mz `numeric` m/z values.
#' @param rules `data.frame` with rules in CAMERA's format (see
#'     [adductRules]).
#'
#' @return `matrix` with one row per m/z value and one column per rule.
#'
#' @export
#'
#' @examples
#' rules <- defaultAdductRules("positive")[1:3, ]
#' impliedMasses(c(170.0445, 192.0275), rules)
impliedMasses <- function(mz, rules) {
    charged <- outer(mz, abs(rules$charge))
    shifted <- sweep(charged, 2, rules$massdiff, "-")
    masses <- sweep(shifted, 2, rules$nmol, "/")
    colnames(masses) <- rules$name
    masses
}

#' Neutral-mass voting within one feature group
#'
#' `neutralMassVote()` is the core of [annotateAdducts()] for a single group
#' of features. Every feature is tried as every quasi-molecular ion
#' (`quasi == 1`). A neutral mass is a hypothesis when at least one other
#' feature is explained by some rule as an ion of the same mass. Duplicate
#' hypotheses (same features with the same rules) are merged, and the rest
#' are ranked by summed rule score (`ips`), then by explained intensity.
#'
#' @param mz `numeric` m/z values of the group's features.
#' @param intensity `numeric` representative intensities of the features.
#' @param rules `data.frame` with rules in CAMERA's format.
#' @param ppm `numeric(1)` relative tolerance in ppm.
#' @param mzabs `numeric(1)` absolute tolerance.
#'
#' @return `data.frame` with one row per hypothesis and member feature and the
#'     columns `feature` (index into `mz`), `hypothesis` (rank), `M`,
#'     `adduct`, `n_ions`, `score` and `intensity`.
#'
#' @export
#'
#' @examples
#' mz <- c(124.0392, 152.0342, 170.0445, 192.0275, 377.0300)
#' neutralMassVote(mz, c(81, 50, 100, 11, 2), defaultAdductRules())
neutralMassVote <- function(mz, intensity, rules, ppm = 10, mzabs = 0.005) {
    if (length(mz) < 2L)
        return(.empty_votes())
    cand <- .nm_candidates(mz, rules)
    anchors <- which(cand$quasi == 1L)
    hyps <- lapply(anchors, .hypothesis_members, cand = cand, ppm = ppm,
                   mzabs = mzabs)
    hyps <- hyps[vapply(hyps, nrow, integer(1)) > 1L]
    if (!length(hyps))
        return(.empty_votes())
    keys <- vapply(hyps, .hypothesis_key, character(1))
    hyps <- hyps[!duplicated(keys)]
    score <- vapply(hyps, .hypothesis_score, numeric(1))
    explained <- vapply(hyps, .explained_intensity, numeric(1),
                        intensity = intensity)
    ord <- order(-score, -explained)
    rows <- mapply(.hypothesis_rows, hyps[ord], seq_along(ord), score[ord],
                   explained[ord], SIMPLIFY = FALSE)
    res <- do.call(rbind, rows)
    rownames(res) <- NULL
    res
}

#' Reduce adduct hypotheses to one per feature
#'
#' `bestHypotheses()` reproduces CAMERA's final choice: within each group it
#' accepts hypotheses from the best ranked down, skipping any hypothesis
#' that contains a feature already explained by an accepted one. Each
#' feature then keeps at most one neutral mass.
#'
#' @param x `data.frame` of hypotheses as returned by [annotateAdducts()] or
#'     [adductHints()].
#'
#' @return The accepted rows of `x`.
#'
#' @export
#'
#' @examples
#' mz <- c(124.0392, 152.0342, 170.0445, 192.0275, 377.0300)
#' hints <- annotateAdducts(mz, AdductVoteParam(),
#'                          intensity = c(81, 50, 100, 11, 2),
#'                          groups = rep("FG.1", 5))
#' bestHypotheses(hints)
bestHypotheses <- function(x) {
    if (!nrow(x))
        return(x)
    x <- x[order(x$group, x$hypothesis), , drop = FALSE]
    key <- paste(x$group, x$hypothesis, sep = "\r")
    keys <- unique(key)
    keep <- logical(length(keys))
    used <- integer()
    ## sequential by nature: each decision depends on the previous ones
    for (k in seq_along(keys)) {
        members <- x$feature[key == keys[k]]
        if (any(members %in% used))
            next
        used <- c(used, members)
        keep[k] <- TRUE
    }
    res <- x[key %in% keys[keep], , drop = FALSE]
    rownames(res) <- NULL
    res
}

## Internal helpers -----------------------------------------------------------

.empty_votes <- function() {
    data.frame(feature = integer(), hypothesis = integer(), M = numeric(),
               adduct = character(), n_ions = integer(), score = numeric(),
               intensity = numeric())
}

.empty_hints <- function() {
    res <- .empty_votes()
    res$group <- character()
    res$rank <- integer()
    res[, .HINT_COLUMNS]
}

.HINT_COLUMNS <- c("feature", "group", "hypothesis", "rank", "M", "adduct",
                   "n_ions", "score", "intensity")

## Long table of every (feature, rule) combination with a positive mass
.nm_candidates <- function(mz, rules) {
    masses <- impliedMasses(mz, rules)
    cand <- data.frame(feature = as.vector(row(masses)),
                       rule = as.vector(col(masses)),
                       M = as.vector(masses))
    cand <- cand[cand$M > 0, , drop = FALSE]
    cand$adduct <- rules$name[cand$rule]
    cand$ips <- rules$ips[cand$rule]
    cand$quasi <- rules$quasi[cand$rule]
    cand
}

## The anchor row of `cand` plus one row per other feature implying the same
## neutral mass (best-scoring rule per feature)
.hypothesis_members <- function(anchor, cand, ppm, mzabs) {
    M0 <- cand$M[anchor]
    tol <- M0 * ppm * 1e-6 + mzabs
    hit <- cand[abs(cand$M - M0) < tol &
                    cand$feature != cand$feature[anchor], , drop = FALSE]
    hit <- hit[order(-hit$ips), , drop = FALSE]
    hit <- hit[!duplicated(hit$feature), , drop = FALSE]
    res <- rbind(cand[anchor, , drop = FALSE], hit)
    res$M0 <- M0
    res
}

.hypothesis_key <- function(h) {
    paste(sort(paste(h$feature, h$adduct)), collapse = "|")
}

.hypothesis_score <- function(h) {
    sum(h$ips)
}

.explained_intensity <- function(h, intensity) {
    sum(intensity[h$feature], na.rm = TRUE)
}

.hypothesis_rows <- function(h, k, score, explained) {
    data.frame(feature = h$feature, hypothesis = k, M = h$M0[1L],
               adduct = h$adduct, n_ions = nrow(h), score = score,
               intensity = explained)
}

## neutralMassVote() for one group, indices mapped back to all features
.vote_group <- function(idx, mz, intensity, rules, ppm, mzabs) {
    res <- neutralMassVote(mz[idx], intensity[idx], rules, ppm, mzabs)
    res$feature <- idx[res$feature]
    res
}

.add_group <- function(x, group) {
    x$group <- rep(group, nrow(x))
    x
}

## Workhorse of the numeric method
.annotate_adducts <- function(mz, intensity, groups, rules, exclude, ppm,
                              mzabs) {
    keep <- which(!exclude)
    by_group <- split(keep, groups[keep])
    votes <- lapply(by_group, .vote_group, mz = mz, intensity = intensity,
                    rules = rules, ppm = ppm, mzabs = mzabs)
    votes <- mapply(.add_group, votes, names(votes), SIMPLIFY = FALSE)
    res <- do.call(rbind, unname(votes))
    if (is.null(res) || !nrow(res))
        return(.empty_hints())
    res$rank <- as.integer(ave(res$hypothesis, res$feature, FUN = rank))
    res <- res[order(res$feature, res$rank), .HINT_COLUMNS, drop = FALSE]
    rownames(res) <- NULL
    res
}
