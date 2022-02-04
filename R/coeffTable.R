#' Table of pairwise relatedness coefficients
#'
#' Creates a data frame containing various relatedness coefficients between all
#' pairs of individuals in a given pedigree.
#'
#' Available coefficients (indicated in `coeffs`) include:
#'
#' * f: The inbreeding coefficient of each pair member. Columns: `f1` and `f2`.
#'
#' * phi: The kinship coefficient. Column: `phi`.
#'
#' * deg: The degree of relationship, as computed by [kin2deg]. Column: `deg`
#'
#' * kappa: The IBD coefficients computed by [kappaIBD]. (These are NA for pairs
#' involving inbred individuals.) Columns: `kappa0`, `kappa1`, `kappa2`.
#'
#' * Delta: The condensed identity coefficients of Jacquard, computed by
#' [condensedIdentity()]. Columns: `D1`, ..., `D9`.
#'
#' @param x A pedigree in the form of a pedtools::ped object.
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param coeffs A character vector containing one or more of the keywords "f",
#'   "phi", "deg", "kappa", "Delta".
#'
#' @return A data frame.
#'
#' @examples
#' # Uncle-nephew pedigree
#' x = addSon(nuclearPed(2), 4)
#'
#' # Complete table
#' coeffTable(x)
#'
#' # Only relevant coefficients
#' coeffTable(x, coeffs = c("phi", "deg", "kappa"))
#'
#' # Only the uncle-nephew pair
#' coeffTable(x, ids = c(3, 6), coeffs = c("phi", "deg", "kappa"))
#' @export
coeffTable = function(x, ids = labels(x), coeffs = c("f", "phi", "deg", "kappa", "Delta")) {
  coeffs = match.arg(coeffs, several.ok = TRUE)

  kappa = kappaIBD(x, ids = ids, simplify = FALSE, inbredAction = 0)
  id1 = kappa$id1
  id2 = kappa$id2

  # Initialise output data frame
  res = data.frame(id1 = id1, id2 = id2)

  if("f" %in% coeffs) {
    inb = inbreeding(x, ids = ids)
    res = cbind(res, f1 = inb[id1], f2 = inb[id2])
  }

  if("phi" %in% coeffs || "deg" %in% coeffs) {
    phiMat = kinship(x)
    phi = phiMat[cbind(id1, id2)]
  }

  if("phi" %in% coeffs) {
    res = cbind(res, phi = phi)
  }

  if("deg" %in% coeffs) {
    res = cbind(res, deg = kin2deg(phi, unrelated = NA_integer_))
  }

  if("kappa" %in% coeffs) {
    res = cbind(res, kappa[, -(1:2)])
  }

  if("Delta" %in% coeffs) {
    delta = condensedIdentity(x, ids = ids, simplify = FALSE, verbose = FALSE)
    res = cbind(res, delta[, -(1:2)])
  }

  res
}




#' Degree of relationship
#'
#' Converts a vector of kinship coefficients to "degrees of relationship", as
#' used by some software for relatedness inference (e.g. KING).
#'
#' The implementation uses the conversion formula \deqn{deg = round(-log2(kin) -
#' 1).}
#'
#' @param kin A vector of kinship coefficients, i.e., numbers in `[0, 1]`.
#' @param unrelated The conversion of unrelated individuals (`kin = 0`).
#'   Mathematically this corresponds to `degree = Inf`, but in some situations
#'   `degree = NA` or something else might be preferable.
#'
#' @return An integer vector of the same length as `kin`.
#' @references KING manual with thresholds for relationship degrees:
#'   <https://www.kingrelatedness.com/manual.shtml>
#'
#' @seealso [kinship()], [coeffTable()]
#' @examples
#' x = cousinPed(1)
#'
#' # Kinship matrix
#' k = kinship(x)
#'
#' # Degrees
#' deg = kin2deg(k)
#' deg
#'
#' # First cousins are 3rd degree
#' stopifnot(deg['7', '8'] == 3)
#'
#' @export
kin2deg = function(kin, unrelated = Inf) {
  deg = round(-log2(kin) - 1)

  if(!identical(unrelated, Inf)) {
    if(length(unrelated) != 1)
      stop2("Argument `unrelated` must have length 1: ", unrelated)
    if(!is.na(unrelated) && is.na(suppressWarnings(as.integer(unrelated))))
      stop2("Argument `unrelated` must be a number or NA: ", unrelated)
    deg[deg == Inf] = unrelated
  }

  deg
}
