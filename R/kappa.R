#' IBD (kappa) coefficients
#'
#' Computes the three IBD coefficients summarising the relationship between two
#' non-inbred individuals.
#'
#' For any pair of non-inbred individuals A and B, their genetic relationship
#' can be summarized by the IBD coefficients \eqn{(\kappa_0, \kappa_1,
#' \kappa_2)}{(\kappa0, \kappa1, \kappa2)}, where \deqn{\kappa_i = P(A and B
#' share i alleles IBD at a random autosomal locus)}
#'
#' The current implementation calls [condensedIdentity()] and returns the three last
#' coefficients in the reverse order.
#'
#' The function checks if any of the `ids` individuals are inbred. If so, a
#' message is printed to the screen, and `c(NA, NA, NA)` is returned.
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object
#' @param ids A character (or coercible to character) of length 2, containing ID
#'   labels of two pedigree members
#' @param ... Arguments passed on to `condensedIdentity`
#'
#' @return The probability vector \eqn{(\kappa_0, \kappa_1, \kappa_2)}{(\kappa0,
#'   \kappa1, \kappa2)}. If any of the two individuals are inbred, `c(NA, NA,
#'   NA)` is returned.
#'
#' @seealso [condensedIdentity()]
#' @examples
#' # Siblings
#' x = nuclearPed(2)
#' k = kappa(x, 3:4)
#' stopifnot(identical(k, c(.25, .5, .25)))
#'
#' # Quad half first cousins
#' x = quadHalfFirstCousins()
#' k = kappa(x, leaves(x))
#' stopifnot(identical(k, c(17/32, 14/32, 1/32)))
#'
#' # Paternal half brothers with 100% inbred father
#' # Genetically indistinguishable from an (outbred) father-son relationship
#' x = halfSibPed()
#' founderInbreeding(x, 1) = 1
#' k = kappa(x, 4:5)
#' stopifnot(identical(k, c(0, 1, 0)))
#'
#' @export
kappa = function(x, ids, ...) {
  j = condensedIdentity(x, ids, ...)
  id1_inbred = any(j[1:4] > .Machine$double.eps)
  id2_inbred = any(j[c(1,2,5,6)] > .Machine$double.eps)
  if(id1_inbred || id2_inbred) {
    message("Inbred individuals: ", toString(ids[c(id1_inbred, id2_inbred)]),
            "\nThe kappa coefficients are undefined for inbred individuals.")
    return(c(NA_real_, NA_real_, NA_real_))
  }
  j[9:7]
}
