furo_hints <- function(nBest = 3) {
    iso <- annotateIsotopes(furo_mz, IsotopologueParam(),
                            intensity = furo_int, groups = furo_grp)
    annotateAdducts(furo_mz, AdductVoteParam(nBest = nBest),
                    intensity = furo_int, groups = furo_grp,
                    exclude = isIsotopologue(iso))
}

test_that("impliedMasses inverts the rules", {
    rules <- defaultAdductRules("positive")
    m <- impliedMasses(170.0445, rules[rules$name == "[M+H]+", ])
    expect_equal(as.numeric(m), furo_M, tolerance = 1e-5)
    m <- impliedMasses(192.0275, rules[rules$name == "[M+Na]+", ])
    expect_equal(as.numeric(m), furo_M, tolerance = 1e-5)
})

test_that("furoylglycine is the best hypothesis", {
    h <- furo_hints()
    expect_identical(colnames(h), .HINT_COLUMNS)
    best <- h[h$hypothesis == 1, ]
    expect_equal(unique(best$M), furo_M, tolerance = 1e-4)
    expect_true(all(c(3L, 4L, 6L) %in% best$feature))
    expect_identical(best$adduct[best$feature == 4L], "[M+H]+")
    expect_identical(best$adduct[best$feature == 6L], "[M+Na]+")
})

test_that("heavy isotopologues are excluded", {
    expect_false(5L %in% furo_hints()$feature)
})

test_that("rank is per feature and starts at 1", {
    h <- furo_hints()
    expect_true(all(tapply(h$rank, h$feature, min) == 1L))
    expect_false(any(duplicated(paste(h$feature, h$rank))))
})

test_that("neutralMassVote handles small inputs", {
    rules <- defaultAdductRules()
    expect_identical(nrow(neutralMassVote(numeric(), numeric(), rules)), 0L)
    expect_identical(nrow(neutralMassVote(170.0445, 1, rules)), 0L)
    ## two unrelated ions
    expect_identical(nrow(neutralMassVote(c(100.0000, 333.3333), c(1, 1),
                                          rules[rules$nmol == 1 &
                                                    rules$charge == 1, ])),
                     0L)
})

test_that("annotateAdducts returns an empty table without hypotheses", {
    h <- annotateAdducts(c(100, 500), AdductVoteParam(), intensity = c(1, 1),
                         groups = c("A", "B"))
    expect_identical(nrow(h), 0L)
    expect_identical(colnames(h), .HINT_COLUMNS)
})

test_that("bestHypotheses keeps one neutral mass per feature", {
    best <- bestHypotheses(furo_hints())
    expect_false(any(duplicated(best$feature)))
    expect_equal(best$M[best$feature == 4L], furo_M, tolerance = 1e-4)
})

test_that("formatAdducts writes the nBest best hints", {
    h <- furo_hints()
    one <- formatAdducts(h, nBest = 1, nFeatures = length(furo_mz))
    expect_length(one, length(furo_mz))
    expect_identical(one[4], "[M+H]+ 169.037")
    expect_identical(one[5], "")
    all_hints <- formatAdducts(h, nBest = Inf, nFeatures = length(furo_mz))
    n_shown <- lengths(strsplit(all_hints[4], " ")) / 2
    expect_identical(as.integer(n_shown),
                     length(unique(sprintf("%.3f", h$M[h$feature == 4L]))))
    expect_error(formatAdducts(h, nBest = 0), "nBest")
})
