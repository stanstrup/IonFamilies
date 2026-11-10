## Generates data/adductRules.rda from the commonMZ package.
##
## commonMZ (https://github.com/stanstrup/commonMZ, MIT licence, same author)
## is not on CRAN or Bioconductor, so its rule tables are shipped here as
## data. Run from the package root after updating commonMZ:
##
##     Rscript inst/scripts/make-adductRules.R
##
## MZ_CAMERA() returns the rules in CAMERA's format, already expanded to one
## to four molecules per ion. warn_clash = FALSE because no rules are
## removed here; users can filter the tables themselves.

library(commonMZ)

adductRules <- list(
    positive = as.data.frame(MZ_CAMERA(mode = "pos", warn_clash = FALSE)),
    negative = as.data.frame(MZ_CAMERA(mode = "neg", warn_clash = FALSE))
)

stopifnot(all(c("name", "nmol", "charge", "massdiff", "quasi", "ips") %in%
                  colnames(adductRules$positive)))

save(adductRules, file = file.path("data", "adductRules.rda"),
     compress = "xz")
