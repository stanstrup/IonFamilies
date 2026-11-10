test_that("empty constructors give valid objects", {
    expect_true(validObject(new("IsotopologueParam")))
    expect_true(validObject(new("AdductVoteParam")))
    expect_s4_class(IsotopologueParam(), "IsotopologueParam")
    expect_s4_class(AdductVoteParam(), "AdductVoteParam")
})

test_that("AdductVoteParam uses the default rules of its polarity", {
    expect_identical(AdductVoteParam()@rules,
                     defaultAdductRules("positive"))
    expect_identical(AdductVoteParam(polarity = "negative")@rules,
                     defaultAdductRules("negative"))
    custom <- defaultAdductRules("positive")[1:5, ]
    expect_identical(AdductVoteParam(rules = custom)@rules, custom)
})

test_that("invalid parameters are rejected", {
    expect_error(IsotopologueParam(ppm = -1), "ppm")
    expect_error(IsotopologueParam(polarity = "both"))
    expect_error(AdductVoteParam(mzabs = c(1, 2)), "mzabs")
    expect_error(AdductVoteParam(nBest = 0), "nBest")
    expect_error(AdductVoteParam(rules = data.frame(name = "[M+H]+")),
                 "lacks")
})

test_that("show methods print", {
    expect_output(show(IsotopologueParam()), "IsotopologueParam")
    expect_output(show(AdductVoteParam(subset = 1:3)), "3 samples")
})

test_that("default rule tables have the required columns", {
    for (pol in c("positive", "negative")) {
        rules <- defaultAdductRules(pol)
        expect_true(all(.RULE_COLUMNS %in% colnames(rules)))
        expect_true(any(rules$quasi == 1))
    }
})
