#' Inbreeding coefficients
#'
#' Compute the inbreeding coefficients of all members of a pedigree. Both
#' autosomal and X-chromosomal coefficients are supported. This function is a
#' simple wrapper of [kinship()]. Note that pedigree founders are allowed to be
#' inbred; see [pedtools::founderInbreeding()] for how to set this up, and see
#' Examples below.
#'
#' The autosomal inbreeding coefficient of a pedigree member is defined as the
#' probability that, at a random autosomal locus, the two alleles carried by the
#' member are identical by descent relative to the pedigree. It follows from the
#' definition that the inbreeding coefficient of a non-founder equals the
#' kinship coefficient of the parents.
#'
#' The implementation here uses [kinship()] to compute the kinship matrix, and
#' computes the inbreeding coefficients from the diagonal, by the formula
#' \deqn{f_a = 2*\phi_{aa} - 1.}{f_a = 2*phi_aa - 1.}
#'
#' The X chromosomal inbreeding coefficient of females are defined (and
#' computed) similarly to the autosomal case above. For males is it always
#' defined as 1.
#'
#' @param x A pedigree in the form of a `ped` object, or a list of such.
#' @param ids A vector of ID labels, or NULL (default).
#' @param Xchrom A logical, indicating if the autosomal (default) or
#'   X-chromosomal inbreeding coefficients should be computed.
#'
#' @return If `ids` has length 1, the inbreeding coefficient of this individual
#'   is returned as a single unnamed number.
#'
#'   Otherwise, the output is a named numeric vector containing the inbreeding
#'   coefficients of the indicated pedigree members (if `ids = NULL`: all).
#'
#' @seealso [kinship()]
#'
#' @examples
#' # Child of half siblings: f = 1/8
#' x = halfCousinPed(0, child = TRUE)
#'
#' # Inbreeding vector
#' inbreeding(x)
#'
#' # Simpler output using the `ids` argument:
#' inbreeding(x, ids = 6)
#'
#'
#' ### X-chromosomal inbreeding ###
#'
#' # Males have inbreeding coefficient 1
#' stopifnot(inbreeding(x, ids = 6, Xchrom = TRUE) == 1)
#'
#' y1 = swapSex(x, ids = 6)  # female child
#' stopifnot(inbreeding(y1, ids = 6, Xchrom = TRUE) == 0)
#'
#' y2 = swapSex(y1, ids = 2) # female ancestor
#' stopifnot(inbreeding(y2, ids = 6, Xchrom = TRUE) == 0.25)
#'
#'
#' ### Inbred founder ###
#'
#' # Mother 100% inbred
#' founderInbreeding(x, ids = 2) = 1
#'
#' inbreeding(x)
#'
#'
#' # Example with selfing and complete inbreeding
#' s = selfingPed(1)
#' founderInbreeding(s, 1) = 1
#' stopifnot(inbreeding(s, ids = 2) == 1)
#'
#' @export
inbreeding = function(x, ids = NULL, Xchrom = FALSE) {

  if(length(ids) == 1) {

    # Quick return if founder
    pars = parents(x, ids)
    if(!length(pars))
      return(founderInbreeding(x, ids = ids, named = FALSE, chromType = if(Xchrom) "x" else "autosomal"))

    # Otherwise: kinship coefficient of parents
    if(Xchrom && getSex(x, ids) == 1L)
      return(1)
    else
      return(kinship(x, ids = pars, Xchrom = Xchrom))
  }

  # Use diagonal of kinship matrix
  kin = kinship(x, Xchrom = Xchrom)
  inb = 2 * diag(kin) - 1

  # X: males are always 1
  if(Xchrom)
    inb[getSex(x) == 1L] = 1

  if(!is.null(ids))
    inb = inb[as.character(ids)]

  inb
}


#' @rdname inbreeding
#' @export
inbreedingX = function(x, ids = NULL) {
  message("This function is deprecated. Use `inbreeding(..., Xchrom = TRUE)` instead.")

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
    inb = inb[as.character(ids)]

  inb
}
