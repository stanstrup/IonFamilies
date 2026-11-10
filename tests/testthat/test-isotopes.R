test_that("isIsotopologue recognises heavy isotopologues only", {
    expect_identical(isIsotopologue(c("[1][M]+", "[1][M+1]+", "[2][M+2]-",
                                      "[M+H]+", "")),
                     c(FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("the 13C isotope of furoylglycine [M+H]+ is found", {
    iso <- annotateIsotopes(furo_mz, IsotopologueParam(),
                            intensity = furo_int, groups = furo_grp)
    expect_length(iso, length(furo_mz))
    expect_identical(iso[4], "[1][M]+")
    expect_identical(iso[5], "[1][M+1]+")
})

test_that("pattern numbers are unique across groups", {
    mz <- c(furo_mz, furo_mz + 100)
    int <- c(furo_int, furo_int)
    grp <- rep(c("A", "B"), each = length(furo_mz))
    iso <- annotateIsotopes(mz, IsotopologueParam(), intensity = int,
                            groups = grp)
    ids <- unique(sub("\\].*", "]", iso[iso != ""]))
    expect_length(ids, 2L)
})

test_that("polarity sets the label suffix", {
    iso <- annotateIsotopes(furo_mz, IsotopologueParam(polarity = "negative"),
                            intensity = furo_int, groups = furo_grp)
    expect_true(all(grepl("-$", iso[iso != ""])))
})

test_that("single features and mismatched lengths are handled", {
    expect_identical(annotateIsotopes(100, IsotopologueParam(),
                                      intensity = 1, groups = "A"), "")
    expect_error(annotateIsotopes(furo_mz, IsotopologueParam(),
                                  intensity = 1, groups = furo_grp),
                 "same length")
})
