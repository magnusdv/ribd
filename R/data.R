#' Basic relationships
#'
#' A data frame containing kinship and kappa coefficients for some commonly used
#' relationships.
#'
#' @format A data frame with 12 rows and 7 columns. The last column (`pos`) is
#'   used internally for placing labels on triangle plots.
#'
"basicRelationships"

#' Jicaque pedigree
#'
#' A data frame describing a pedigree from the Jicaque tribe, studied by Chapman
#' and Jacquard (1971).
#'
#' @format A data frame with 22 rows and four columns:
#'
#'  * `id`  : individual ID
#'  * `fid` : father's ID (or 0 if not included)
#'  * `mid` : mother's ID (or 0 if not included)
#'  * `sex` : Gender codes, where 1 = male and 2 = female
#'
#' @references Chapman, A.M and Jacquard, A. (1971). Un isolat d'Amerique
#'   Centrale: les Indiens Jicaques de Honduras. In Genetique et Population.
#'   Paris: Presses Universitaires de France.
"jicaque"
