#' Two-locus inbreeding
#'
#' Computes the two-locus inbreeding coefficient of a pedigree member, for a
#' given recombination rate.
#'
#' Let A be a pedigree member, and L1, L2 two autosomal loci with recombination
#' rate \eqn{\rho}. The two-locus inbreeding coefficient \eqn{f_{11}(\rho)} is
#' defined as the probability that A is autozygous at both L1 and L2
#' simultaneously.
#'
#' As in the one-locus case, the two-locus inbreeding coefficient of A equals
#' the two-locus kinship coefficient of the parents.
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object.
#' @param id The ID label of a pedigree member.
#' @param rho A numeric vector of recombination rates; all entries must be in
#'   the interval \eqn{[0, 0.5]}.
#' @param verbose A logical.
#' @param debug A logical. If TRUE, detailed messages are printed during the
#'   recursion process.
#'
#' @seealso [twoLocusKinship()], [twoLocusIBD()], [twoLocusIdentity()]
#'
#' @references Weir & Cockerham (1969). _Pedigree mating with two linked loci_.
#'   Genetics, 61:923-940.
#'
#'
#' @examples
#'
#' ###################################################
#' # Reproducing an example of Weir & Cockerham (1969)
#' ###################################################
#'
#' # Pedigree
#' x = nuclearPed(2, sex = 1:2) |>
#'   addDaughter(3:4) |>
#'   addSon(c(3,5)) |>
#'   addDaughter(5:6) |>
#'   relabel(new = strsplit("GHDECBA","")[[1]])
#'
#' plot(x)
#'
#' # The two-locus inbreeding of A
#' twoLocusPlot(list(ped = x, ids = "A"), coeff = "inb")
#'
#' # W&C formula (expressed by linkage parameter a = 1-2*rho)
#' rho = seq(0, 0.5, length = 11)
#' a = 1 - 2*rho
#' WC = (128 + 10*a + 36*a^2 + 47*a^3 + 20*a^4 + 10*a^5 + 4*a^6 + a^7)/512
#'
#' points(rho, WC, col = 2)
#'
#' @export
twoLocusInbreeding = function(x, id, rho, verbose = FALSE, debug = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(id) != 1) stop2("Argument `id1` must have length 1")
  if(!is.numeric(rho))
    stop2("Argument `rho` must be numeric")
  if(any(rho < 0 | rho > 0.5))
    stop2("Argument `rho` cannot have entries outside the interval [0, 0.5]: ", rho[rho < 0 | rho > 0.5])
  if(debug && length(rho) > 1)
    stop2("Debugging mode is only allowed when `rho` has length 1")

  pars = parents(x, id)

  # Quick return if founder
  if(!length(pars)) {
    finb = founderInbreeding(x, ids = id, named = FALSE, chromType = "autosomal")
    if(finb > 0 && finb < 1)
      stop2("Two-locus inbreeding is not well-defined under incomplete founder inbreeding")
    return(finb)
  }

  twoLocusKinship(x, ids = pars, rho = rho, verbose = verbose, debug = debug)
}
