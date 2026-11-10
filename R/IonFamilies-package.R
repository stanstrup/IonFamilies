#' IonFamilies: isotope and adduct annotation of LC-MS feature groups
#'
#' One compound gives several LC-MS features: its quasi-molecular ion, other
#' adducts, in-source fragments, and the isotopes of all of these. This
#' package annotates that *ion family* within groups of features, typically
#' the groups made with [MsFeatures::groupFeatures()] on an
#' [xcms::XcmsExperiment].
#'
#' The workflow has two steps:
#'
#' 1. [annotateIsotopes()] with an [IsotopologueParam()] finds isotope
#'    patterns with [MetaboCoreUtils::isotopologues()] and labels them in
#'    CAMERA syntax (`[3][M]+`, `[3][M+1]+`).
#' 2. [annotateAdducts()] with an [AdductVoteParam()] finds adducts and
#'    fragments by neutral-mass voting with a rule list
#'    ([defaultAdductRules()]). All supported hypotheses are kept and ranked
#'    ([adductHints()]); the best `nBest` are written as text
#'    ([formatAdducts()]).
#'
#' Annotations are *hints*: they suggest what to check, they do not identify
#' a compound.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom methods new is validObject setClass setGeneric setMethod
#'     prototype setValidity show slotNames
#' @importFrom stats median ave
#' @importFrom MetaboCoreUtils isotopologues
#' @importFrom MsFeatures featureGroups
#' @importFrom xcms featureDefinitions featureDefinitions<- featureValues
#' @importFrom S4Vectors metadata metadata<-
NULL
