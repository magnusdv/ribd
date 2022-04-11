#' Karigl's generalised kinship coefficients
#'
#' Compute generalised kinship coefficients, as defined by Karigl (1981),
#' involving up to 4 pedigree members. The founders may be inbred; see Examples.
#'
#' The function `generalisedKinship3()` computes the generalised kinship
#' coefficient of three (not necessarily distinct) members `a`, `b` and `c`,
#' defined as the probability that if a random allele is chosen from each of
#' them, they are all identical by descent.
#'
#' The function `generalisedKinship4()` computes the generalised kinship
#' coefficient of four individuals, defined similarly to the above.
#'
#' The function `generalisedKinship22()` computes the generalised kinship
#' coefficient of two pairs of members, defined as the probability that in both
#' pairs simultaneously, random alleles chosen from the two individuals are IBD.
#'
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#' @param ids A vector of ID labels, of length 3 for `generalisedKinship3()` and
#'   4 for `generalisedKinship4()` and `generalisedKinship22()`.
#' @param sparse A positive integer, indicating the pedigree size limit for
#'   using sparse arrays. If NA, a default limit of 50 is used.
#' @param chromType Either "autosomal" or "x".
#' @param verbose A logical.
#'
#' @return A numeric of length 1.
#'
#' @seealso [kinship()], [condensedIdentity()],
#'   [condensedIdentityX()]
#' @examples
#' # Generalised kinship between three siblings
#' x = nuclearPed(3)
#' phi3 = generalisedKinship3(x, ids = 3:5)
#'
#' # Recalculate if the father is 100% inbred
#' founderInbreeding(x, 1) = 1
#' phi3_inbred = generalisedKinship3(x, ids = 3:5)
#'
#' stopifnot(phi3 == 1/16, phi3_inbred == 1/8 + 1/32)
#'
#' @name generalised_karigl
NULL

#' @rdname generalised_karigl
#' @export
generalisedKinship3 = function(x, ids, sparse = NA, chromType = "autosomal", verbose = FALSE) {

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = chromType, verbose = verbose)

  # Compute
  res = phi3(ids_int[1], ids_int[2], ids_int[3], chromType, mem)

  if(verbose)
    printCounts(mem)

  res
}

#' @rdname generalised_karigl
#' @export
generalisedKinship4 = function(x, ids, sparse = NA, chromType = "autosomal", verbose = FALSE) {

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = chromType, verbose = verbose)

  res = phi4(ids_int[1], ids_int[2], ids_int[3], ids_int[4], chromType, mem)
  if(verbose)
    printCounts(mem)

  res
}

#' @rdname generalised_karigl
#' @export
generalisedKinship22 = function(x, ids, sparse = NA, chromType = "autosomal", verbose = FALSE) {
  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = chromType, verbose = verbose)

  res = phi22(ids_int[1], ids_int[2], ids_int[3], ids_int[4], chromType, mem)
  if(verbose)
    printCounts(mem)

  res
}

