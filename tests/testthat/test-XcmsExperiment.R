test_that("feature groups are required", {
    skip_if_not_installed("faahKO")
    x <- xcms::loadXcmsData("xmse")
    expect_error(annotateIsotopes(x, IsotopologueParam()), "groupFeatures")
    expect_error(annotateAdducts(x, AdductVoteParam()), "groupFeatures")
})

test_that("annotation is stored in featureDefinitions and metadata", {
    x <- grouped_xmse()
    x <- annotateIsotopes(x, IsotopologueParam())
    x <- annotateAdducts(x, AdductVoteParam(nBest = 2))
    fd <- xcms::featureDefinitions(x)
    expect_true(all(c("isotopes", "adduct", "adduct_hints", "pcgroup") %in%
                        colnames(fd)))
    expect_type(fd$isotopes, "character")
    expect_type(fd$adduct, "character")
    expect_type(fd$pcgroup, "integer")
    expect_length(fd$adduct_hints, nrow(fd))
    params <- S4Vectors::metadata(x)$IonFamilies
    expect_s4_class(params$isotopes, "IsotopologueParam")
    expect_s4_class(params$adducts, "AdductVoteParam")
})

test_that("pcgroup numbers follow the feature groups", {
    x <- annotateAdducts(grouped_xmse(), AdductVoteParam())
    fd <- xcms::featureDefinitions(x)
    expect_identical(length(unique(fd$pcgroup)),
                     length(unique(fd$feature_group)))
    expect_true(all(tapply(fd$pcgroup, fd$feature_group, n_unique) == 1L))
})

test_that("adductHints round-trips the stored hypotheses", {
    x <- grouped_xmse()
    x <- annotateIsotopes(x, IsotopologueParam())
    x <- annotateAdducts(x, AdductVoteParam())
    h <- adductHints(x)
    expect_identical(colnames(h), .HINT_COLUMNS)
    fd <- xcms::featureDefinitions(x)
    expect_identical(sort(unique(h$feature)),
                     which(vapply(fd$adduct_hints, nrow, integer(1)) > 0L))
    expect_false(any(isIsotopologue(fd$isotopes[h$feature])))
})

test_that("formatAdducts changes nBest without recomputing", {
    x <- annotateAdducts(grouped_xmse(), AdductVoteParam(nBest = 1))
    x2 <- formatAdducts(x, nBest = Inf)
    expect_identical(adductHints(x), adductHints(x2))
    expect_identical(S4Vectors::metadata(x2)$IonFamilies$adducts@nBest, Inf)
    n1 <- lengths(strsplit(xcms::featureDefinitions(x)$adduct, " "))
    n2 <- lengths(strsplit(xcms::featureDefinitions(x2)$adduct, " "))
    expect_true(all(n2 >= n1))
    expect_true(any(n2 > n1))
})

test_that("annotation survives filterFeatureDefinitions", {
    x <- annotateAdducts(grouped_xmse(), AdductVoteParam())
    y <- xcms::filterFeatureDefinitions(x, 1:20)
    fd <- xcms::featureDefinitions(y)
    expect_identical(nrow(fd), 20L)
    expect_true(all(c("adduct", "adduct_hints", "pcgroup") %in%
                        colnames(fd)))
    expect_s4_class(S4Vectors::metadata(y)$IonFamilies$adducts,
                    "AdductVoteParam")
})

test_that("adductHints needs an adduct annotation", {
    expect_error(adductHints(grouped_xmse()), "annotateAdducts")
})
