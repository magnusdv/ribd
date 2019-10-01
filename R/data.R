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


#' Multi-person IBD patterns
#'
#' A list of length 6. Each list entry is a matrix of IBD patterns with some
#' extra information in additional columns. The N'th entry contains all
#' condensed patterns of N individuals.
#'
#' @format A list of 6 matrices. The numbers of rows are respectively 1, 3, 16, 139,
#'   1750 and 29388.
"MULTIPATTERNS_NONINBRED"
