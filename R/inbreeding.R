#' Inbreeding coefficients
#'
#' Compute the inbreeding coefficients of all members of a pedigree. These are
#' simple wrappers of [kinship()] and [kinshipX()]. The founders may be inbred;
#' see [pedtools::founderInbreeding()] for how to set this up.
#'
#' The autosomal inbreeding coefficient of a pedigree member is defined as the
#' probability that, at a random autosomal locus, the two alleles carried by the
#' member are identical by descent relative to the pedigree. It follows from the
#' definition that the inbreeding coefficient of a member equals the kinship
#' coefficient of the parents.
#'
#' The X chromosomal inbreeding coefficient of an female member is defined
#' similarly to the autosomal case above. For males is it always 1.
#'
#' The inbreeding coefficients are computed from the diagonal of the kinship
#' matrix, by the formula \deqn{f_a = 2*\phi_{aa} - 1.}{f_a = 2*phi_aa - 1.}
#'
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#' @param ids A vector of ID labels, or NULL (default).
#' @return If `ids` has length 1, the inbreeding coefficient of this individual
#'   is returned as a single unnamed number.
#'
#'   Otherwise, the output is a named numeric vector containing the inbreeding
#'   coefficients of the indicated pedigree members (if `ids = NULL`: all).
#'
#' @seealso [kinship()]
#' @examples
#' # Child of half siblings: f = 1/8
#' x = halfCousinPed(0, child = TRUE)
#' inbreeding(x)
#'
#' # If the father is 100% inbred, the inbreeding coeff of the child doubles
#' fa = commonAncestors(x, 4:5) # robust to label change
#' founderInbreeding(x, fa) = 1
#'
#' inbreeding(x)
#'
#' # Simpler output using the `ids` argument:
#' inbreeding(x, ids = 6)
#'
#' ### X-chromosomal inbreeding coefficients ###
#' # These depend on the genders in the pedigree.
#' # To exemplify, we consider a child of half siblings.
#'
#' xPat = halfSibPed(sex2 = 2) # paternal half sibs
#' xPat = addChildren(xPat, father = 4, mother = 5, nch = 1, sex = 2)
#' stopifnot(inbreedingX(xPat, ids = 6) == 0)
#'
#' # Change to maternal half sibs => coeff becomes 1/4.
#' xMat = swapSex(xPat, 1)
#' stopifnot(inbreedingX(xMat, ids = 6) == 0.25)
#'
#' # Example with selfing and complete inbreeding
#' s = selfingPed(1)
#' founderInbreeding(s, 1) = 1
#' inbreeding(s, ids = 2)
#'
#' @export
inbreeding = function(x, ids = NULL) {

  if(length(ids) == 1) {

    # Quick return if founder
    pars = parents(x, ids)
    if(!length(pars))
      return(founderInbreeding(x, ids = ids, named = FALSE))

    # Otherwise: kinship coefficient of parents
    return(kinship(x, ids = pars))
  }

  # Use diagonal of kinship matrix
  kin = kinship(x)
  inb = 2 * diag(kin) - 1

  if(!is.null(ids))
    inb = inb[ids]

  inb
}


#' @rdname inbreeding
#' @export
inbreedingX = function(x, ids = NULL) {

  if(length(ids) == 1) {

    # Quick return if founder
    pars = parents(x, ids)
    if(!length(pars))
      return(founderInbreeding(x, ids = ids, named = FALSE, chromType = "X"))

    # Otherwise: kinship coefficient of parents
    return(kinshipX(x, ids = pars))
  }

  # Use diagonal of kinship matrix
  kin = kinshipX(x)
  inb = 2 * diag(kin) - 1

  if(!is.null(ids))
    inb = inb[ids]

  inb
}
