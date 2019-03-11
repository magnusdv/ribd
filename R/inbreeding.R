#' Inbreeding coefficients
#'
#' Compute the inbreeding coefficients of all members of a pedigree. These are
#' wrappers of [kinship()] and [kinshipX()] which do the main work. The founders
#' may be inbred; see [pedtools::founderInbreeding()] for how to set this up.
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
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#'
#' @return A vector of length `pedsize(x)`.
#'
#' @seealso [kinship()]
#' @examples
#' # Child of half siblings: f = 1/8
#' x = halfCousinPed(0, child = TRUE)
#' inbreeding(x)
#'
#' # If the father is 100% inbred, the inbreeding coeff of the child doubles
#' founderInbreeding(x, 1) = 1
#' inbreeding(x)
#'
#' # The X inbreeding coefficients depend on the gender distribution
#' x.pat = halfCousinPed(0, child = TRUE) # child of paternal half sibs
#' inbreedingX(x.pat) # all zero
#'
#' x.mat = swapSex(x.pat, 1) # change to maternal half sibs
#' inbreedingX(x.mat) # child now has f_X = 1/4
#'
#' @export
inbreeding = function(x) {
  kin.matrix = kinship(x)

  # Initialize result vector
  inb = numeric(pedsize(x))
  names(inb) = labels(x)

  # Fill in founder inbreeding
  inb[founders(x)] = founderInbreeding(x)

  # Non-founders
  inb[nonfounders(x)] = kin.matrix[cbind(x$FIDX, x$MIDX)] # founders -> kin[0,0] -> excluded

  inb
}

#' @rdname inbreeding
#' @export
inbreedingX = function(x) {
  kin.matrix = kinshipX(x)

  # Initialize result vector
  inb = numeric(pedsize(x))
  names(inb) = labels(x)

  # Fill in founder inbreeding
  # inb[founders(x)] = founderInbreeding(x)

  # Non-founders
  inb[nonfounders(x)] = kin.matrix[cbind(x$FIDX, x$MIDX)] # founders -> kin[0,0] -> excluded

  inb
}
