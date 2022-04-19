#' Table of pairwise relatedness coefficients
#'
#' Creates a data frame containing various relatedness coefficients between all
#' pairs of individuals in a given pedigree.
#'
#' Available coefficients (indicated in `coeff`) include:
#'
#' * `f`: The inbreeding coefficient of each pair member. Columns: `f1` and
#' `f2`.
#'
#' * `phi`: The kinship coefficient. Column: `phi`.
#'
#' * `deg`: The degree of relationship, as computed by [kin2deg]. Column: `deg`
#'
#' * `kappa`: The IBD coefficients computed by [kappaIBD]. (These are NA for
#' pairs involving inbred individuals.) Columns: `k0`, `k1`, `k2`.
#'
#' * `identity`: The 9 condensed identity coefficients of Jacquard, computed by
#' [identityCoefs()]. Columns: `D1`, ..., `D9`.
#'
#' * `detailed`: The detailed identity coefficients of Jacquard, computed by
#' `identityCoefs(..., detailed = TRUE)`. Columns: `d1`, ..., `d15`.
#'
#' @param x A pedigree in the form of a pedtools::ped object.
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param coeff A character vector containing one or more of the keywords "f",
#'   "phi", "deg", "kappa", "identity", "detailed".
#' @param self A logical indicating if self-relationships should be included.
#'   Default: FALSE.
#' @param Xchrom A logical indicating if the coefficients should be autosomal
#'   (default) or X-chromosomal. If `Xchrom = NA`, both sets are included.
#'
#' @return A data frame with one row for each pair of individuals. The first two
#'   columns are characters named `id1` and `id2`, while remaining columns are
#'   numeric. For columns containing X-chromosomal coefficients, their names are
#'   prefixed with "X-".
#'
#'   If "f" (inbreeding) is the only coefficient, the data frame has one row per
#'   individual, and the first column is named `id`.
#'
#'   Note: If `x` has members with unknown members, all X-chromosomal
#'   coefficients are NA.
#'
#'
#' @examples
#' # Uncle-nephew pedigree
#' x = addSon(nuclearPed(2), 4)
#'
#' # Complete table
#' coeffTable(x)
#'
#' # Only relevant coefficients
#' coeffTable(x, coeff = c("phi", "deg", "kappa"))
#'
#' # Only the uncle-nephew pair
#' coeffTable(x, ids = c(3, 6), coeff = c("phi", "deg", "kappa"))
#'
#' # X-chromosomal coefficients
#' coeffTable(x, Xchrom = TRUE)
#'
#' # Both autosomal and X
#' coeffTable(x, Xchrom = NA)
#'
#' @export
coeffTable = function(x, ids = labels(x), coeff = c("f", "phi", "deg", "kappa", "identity", "detailed"),
                      self = FALSE, Xchrom = FALSE) {

  if(!is.ped(x)) {
    if(is.pedList(x))
      stop2("Lists of pedigrees are currently not supported in `coeffTable()`")
    else
      stop2("Illegal input. Expected `ped`, received: ", class(x))
  }

  if(anyDuplicated(ids))
    stop2("Duplicated ID label: ", unique.default(ids[duplicated(ids)]))

  coeff = match.arg(coeff, several.ok = TRUE)

  # If only inbreeding, treat separately and return
  if(length(coeff) == 1 && coeff == "f"){
    res = data.frame(id = ids)
    if(is.na(Xchrom) || !Xchrom)
      res = cbind(res, f = inbreeding(x, ids))
    if(is.na(Xchrom) || Xchrom) {
      xf = if(any(x$SEX == 0)) rep(NA_real_, length(ids)) else inbreeding(x, ids, Xchrom = TRUE)
      res = cbind(res, `X-f` = xf)
    }
    return(res)
  }

  # Xchrom = NA: both autosomal and X
  if(is.na(Xchrom)) {
    aut = coeffTable(x, ids, coeff = coeff, Xchrom = FALSE, self = self)
    xchr = coeffTable(x, ids, coeff = coeff, Xchrom = TRUE, self = self)
    return(cbind(aut, xchr[, -(1:2)]))
  }


  # If Xchrom & unknown sex: Fill with NAs
  if(Xchrom && any(x$SEX == 0)) {
    # horrible hack: compute autosomal version to get proper names and dimensions
    res = coeffTable(x, ids, coeff = coeff, self = self)
    res[, -(1:2)] = NA_real_
    names(res)[-(1:2)] = paste0("X-", names(res)[-(1:2)])
    return(res)
  }

  # Initialise output data frame
  idmat = .idPairs(ids, self = self, as = "character", returnList = FALSE)
  res = data.frame(id1 = idmat[,1], id2 = idmat[,2])

  if("f" %in% coeff) {
    inb = inbreeding(x, ids = ids, Xchrom = Xchrom)
    res = cbind(res, f1 = inb[res$id1], f2 = inb[res$id2])
  }

  # Kinship matrix: needed if phi and/or deg are included
  if("phi" %in% coeff || "deg" %in% coeff) {
    phiMat = kinship(x, Xchrom = Xchrom)
    phi = phiMat[cbind(res$id1, res$id2)]
  }

  if("phi" %in% coeff) {
    res = cbind(res, phi = phi)
  }

  if("deg" %in% coeff) {
    res = cbind(res, deg = kin2deg(phi, unrelated = NA_integer_))
  }

  if("kappa" %in% coeff) {
    kappa = kappaIBD(x, ids = ids, simplify = FALSE, inbredAction = 0, Xchrom = Xchrom)
    names(kappa)[3:5] = paste0("k", 0:2) # abbreviate prefix from "kappa" to "k"
    res = merge(res, kappa, by = 1:2, all.x = TRUE)
  }

  if("identity" %in% coeff) {
    delta = identityCoefs(x, ids = ids, detailed = FALSE, simplify = FALSE, self = self,
                          Xchrom = Xchrom, verbose = FALSE)
    res = cbind(res, delta[, -(1:2)])
  }

  if("detailed" %in% coeff) {
    delta = identityCoefs(x, ids = ids, detailed = TRUE, simplify = FALSE, self = self,
                          Xchrom = Xchrom, verbose = FALSE)
    res = cbind(res, delta[, -(1:2)])
  }

  if(Xchrom)
    names(res)[-(1:2)] = paste0("X-", names(res)[-(1:2)])

  res
}




#' Degree of relationship
#'
#' Converts a vector of kinship coefficients to "degrees of relationship", as
#' used by some software for relatedness inference (e.g. KING).
#'
#' The implementation uses the conversion formula \deqn{deg = round(-log2(kin) -
#' 1).}
#' The first degrees correspond to the following approximate kinship ranges:
#'
#' * `[0.354, 1]`: 0th degree (MZ twins or duplicates)
#'
#' * `[0.177, 0.354)`: 1st degree (parent-offspring, full siblings)
#'
#' * `[0.0884, 0.177)`: 2nd degree (half sibs, grandparent-grandchild, avuncular)
#'
#' * `[0.0442, 0.0884)` 3rd degree (half-avuncular, first cousins, great-grandparent etc)
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

  # Push negative degrees (inbred MZ) up to 0
  deg[deg < 0] = 0

  # Modify entries for unrelated?
  if(!identical(unrelated, Inf)) {
    if(length(unrelated) != 1)
      stop2("Argument `unrelated` must have length 1: ", unrelated)
    if(!is.na(unrelated) && is.na(suppressWarnings(as.integer(unrelated))))
      stop2("Argument `unrelated` must be a number or NA: ", unrelated)
    deg[deg == Inf] = unrelated
  }

  deg
}
