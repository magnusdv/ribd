#' Two-locus IBD coefficients (EXPERIMENTAL)
#'
#' Computes the two-locus IBD coefficient \eqn{\kappa_{1,1}}{k_11} of a pair of
#' pedigree members, at a given recombination rate.
#'
#' Let A, B be two pedigree members, and L1, L2 two loci with a given
#' recombination rate r. The two-locus IBD coefficient
#' \eqn{\kappa_{1,1}(\rho)}{k_11(\rho)} is defined as the probability that A and B
#' have IBD status 1 at both L1 and L2 simultaneously. Note that it is not
#' required that the IBD alleles are in cis (or in trans for that matter), in
#' either individual.
#'
#' The method of computation is via the two-locus kinship coefficients involving
#' the parents of A and B.
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object.
#' @param ids A character (or coercible to character) containing ID labels of
#'   two pedigree members.
#' @param rho A number in the interval \eqn{[0, 0.5]}; the recombination rate
#'   between the two loci.
#'
#' @return A single numeric.
#'
#' @seealso [twoLocusKinship()]
#'
#' @examples
#'
#' ######################################
#' # A classic example of three relationships with the same one-locus IBD
#' # coefficients, but different two-locus IBD coefficients. As a consequence,
#' # these relationships cannot be separated using unlinked markers, but they
#' # are separable (in theory) using pairs of linked markers.
#' ######################################
#' peds = list(
#'     GrandParent  = list(ped = linearPed(2),    ids = c(1,5)),
#'     HalfSib      = list(ped = halfSibPed(),    ids = c(4,5)),
#'     Uncle        = list(ped = cousinPed(0, 1), ids = c(3,6)))
#'
#' # Recombination values
#' rseq = seq(0, 0.5, length = 11)
#'
#' # Compute two-locus IBD coefficients (k11)
#' kvals = sapply(peds, function(x)
#'   sapply(rseq, function(r) twoLocusIBD(x$ped, x$ids, r)))
#'
#' # Plot
#' matplot(rseq, kvals, type = "l", xlab = "Recombination rate", ylab = "",
#'         main = expression(paste("Two-locus IBD:  ", kappa[`1,1`])))
#'
#' legend("topright", names(peds), col = 1:3, lty = 1:3)
#'
#' #########################
#' # A more complex pedigree
#' #########################
#' y = addSon(cousinPed(0, child = TRUE), 5)
#' ids = c(1, 7) # neither is inbred
#'
#' # Exact two-locus k11
#' k11 = sapply(rseq, function(r) twoLocusIBD(y, ids, r))
#' plot(rseq, k11, type="l", ylim = c(0, 0.5))
#'
#' \dontrun{
#' # Check by simulation (increase Nsim!)
#' library(ibdsim2)
#' k11.sim = sapply(rseq, function(r)
#'                  estimateTwoLocusIBD(y, ids, r, Nsim = 100)['ibd1', 'ibd1'])
#' points(rseq, k11.sim, col = 2)
#' }
#'
#' @export
twoLocusIBD = function(x, ids, rho) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must have length exactly 2")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  ### Sort ids so that ids[1] comes first in the internal ordering
  # if(ids_int[1] > ids_int[2]) {
  #  ids[] = ids[2:1]
  #  ids_int[] = ids_int[2:1]
  # }
  # directDescend = ids_int[1] %in% ancestors(x, ids_int[2], internal = T)

  # Check that none of the ids are inbred
  if(any(inbreeding(x)[ids_int] > 0))
    stop2("IBD coefficients are not defined for pairs of inbred individuals")

  # One-locus IBD coefficients
  kap = kappaIbd(x, ids)

  # If kappa2 = 0
  if(kap[3] == 0) {
    p11.hh = twoLocusKinship(x, ids, rho, recombinants = c(F,F))
    p11.rh = twoLocusKinship(x, ids, rho, recombinants = c(T,F))
    p11.hr = twoLocusKinship(x, ids, rho, recombinants = c(F,T))
    p11.rr = twoLocusKinship(x, ids, rho, recombinants = c(T,T))

    # Compute k11 as sum of four probs: cc = cis/cis; ct = cis/trans etc
    # If rho = 0, then only the cc part is nonzero
    k11.cc = p11.hh * 4/(1-rho)^2

    if(rho == 0)
      return(k11.cc)

    k11.ct = p11.hr * 4/(rho*(1-rho))
    k11.tc = p11.rh * 4/(rho*(1-rho))
    k11.tt = p11.rr * 4/(rho^2)

    return(k11.cc + k11.ct + k11.tc + k11.tt)
  }
  stop2("Relationships with kappa2 > 0 are not implemented yet")

  id1 = ids[1]
  id2 = ids[2]
  fa1 = father(x, id1)
  fa2 = father(x, id2)
  mo1 = mother(x, id1)
  mo2 = mother(x, id2)

  phi = kinship(x)

  # Matrix of two-locus kinship for all pairs
  # TODO: only include ids up to max(internalID(ids))
  phi11 = matrix(NA, ncol = pedsize(x), nrow = pedsize(x),
                 dimnames = list(labels(x), labels(x)))

  phi11.df = twoLocusKinship(x, ids=labels(x), rho)
  int1 = internalID(x, phi11.df$id1)
  int2 = internalID(x, phi11.df$id2)
  phi11[cbind(int1, int2)] = phi11[cbind(int2, int1)] = phi11.df$phi2

  # Derive similar matrices of phi10, phi01 and phi00
  phi01 = phi10 = phi - phi11
  phi00 = 1 - 2*phi + phi11

  # The following probabilities are known:
  k22.h = phi11[fa1, fa2] * phi11[mo1, mo2] + phi11[fa1, mo2] * phi11[mo1, fa2]
  k21.h = phi11[fa1, fa2] * phi10[mo1, mo2] + phi11[fa1, mo2] * phi10[mo1, fa2] +
          phi11[mo1, fa2] * phi10[fa1, mo2] + phi11[mo1, mo2] * phi10[fa1, fa2]
  k11.cc = phi11[fa1, fa2] * phi00[mo1, mo2] + phi11[fa1, mo2] * phi00[mo1, fa2] +
           phi11[mo1, fa2] * phi00[fa1, mo2] + phi11[mo1, mo2] * phi00[fa1, fa2]

  phi11.rr = twoLocusKinship(x, ids=ids, rho, recombinants = c(T,T))
  k11.tt = phi11.rr * 4/rho^2 - 2*k22.h - 2*k21.h

  return(c(k22.h = k22.h, k21.h = k21.h, k11.cc = k11.cc, k11.tt = k11.tt))
}
