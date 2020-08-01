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
#' @param id Either a single ID label, or NULL (default).
#' @return If `id` is NULL, the output is a named numeric vector of length
#'   `pedsize(x)`, containing the inbreeding coefficients of each pedigree
#'   member.  is returned.
#'
#'   If `id` is the label of a pedigree member, the inbreeding coefficient of
#'   this individual is returned unnamed.
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
#' # Simpler output using the `id` argument:
#' inbreeding(x, id = 6)
#'
#' ### X-chromosomal inbreeding coefficients ###
#' # These depend on the genders in the pedigree.
#' # To exemplify, we consider a child of half siblings.
#'
#' xPat = halfSibPed(sex2 = 2) # paternal half sibs
#' xPat = addChildren(xPat, father = 4, mother = 5, nch = 1, sex = 2)
#' stopifnot(inbreedingX(xPat, id = 6) == 0)
#'
#' # Change to maternal half sibs => coeff becomes 1/4.
#' xMat = swapSex(xPat, 1)
#' stopifnot(inbreedingX(xMat, id = 6) == 0.25)
#'
#' # Example with selfing and complete inbreeding
#' s = selfingPed(1)
#' founderInbreeding(s, 1) = 1
#' inbreeding(s, id = 2)
#'
#' @export
inbreeding = function(x, id = NULL) {

  if(!is.null(id)) {
    if(length(id) != 1)
      stop2("When `id` is not NULL, it must be a single ID label")

    # Quick return if founder
    pars = parents(x, id)
    if(!length(pars))
      return(founderInbreeding(x, ids = id, named = FALSE))

    # Otherwise: kinship coefficient of parents
    return(kinship(x, ids = pars))
  }

  # If `id = NULL`: use diagonal of kinship matrix
  kin = kinship(x)
  2 * diag(kin) - 1
}


#' @rdname inbreeding
#' @export
inbreedingX = function(x, id = NULL) {

  if(!is.null(id)) {
    if(length(id) != 1)
      stop2("When `id` is not NULL, it must be a single ID label")

    # Quick return if founder
    pars = parents(x, id)
    if(!length(pars))
      return(founderInbreeding(x, ids = id, named = FALSE, chromType = "X"))

    # Otherwise: kinship coefficient of parents
    return(kinshipX(x, ids = pars))
  }

  # If `id = NULL`: use diagonal of kinship matrix
  kin = kinshipX(x)
  2 * diag(kin) - 1
}
