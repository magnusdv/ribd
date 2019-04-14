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
#' @return
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
#'     GrandP  = list(ped = linearPed(2),    ids = c(1,5)))
#'   # Uncle   = list(ped = cousinPed(0, 1), ids = c(3,6)),  # Not yet
#'   # HalfSib = list(ped = halfSibPed(),    ids = c(4,5)))  # Not yet
#'
#' # Recombination values
#' rseq = seq(0, 0.5, length = 11)
#'
#' # Compute two-locus IBD coefficients (k11)
#' kvals = sapply(peds, function(x)
#'   sapply(rseq, function(r) twoLocusIBD(x$ped, x$ids, r)))
#'
#' # Plot
#' matplot(rseq, kvals, type = "l")
#' legend("topright", names(peds), col = 1:3, lty = 1:3)
#'
#' #########################
#' # A more complex pedigree
#' #########################
#' y = addSon(cousinPed(0, child = TRUE), 5)
#' ids = c(1, 7) # neither is inbred
#'
#' # TODO: Check if the following is correct
#' sapply(rseq, function(r) twoLocusIBD(y, ids, r))
#'
#' @export
twoLocusIBD = function(x, ids, rho) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must have length exactly 2")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  # Sort ids so that ids[1] comes first in the internal ordering
  ids_int = internalID(x, ids)
  if(ids_int[1] > ids_int[2]) {
    ids[] = ids[2:1]
    ids_int[] = ids_int[2:1]
  }

  # Check that none of the ids are inbred
  if(any(inbreeding(x)[ids_int] > 0))
    stop2("IBD coefficients are not defined for pairs of inbred individuals")

  # Is ids[2] a direct descendant of ids[1]?
  directDescend = ids_int[1] %in% ancestors(x, ids_int[2], internal = T)

  # In descendant case, use special formula
  if(directDescend) {
    p11_nr_nr = twoLocusKinship(x, ids, rho, recombinants = c(F,F))
    p11_r_nr = twoLocusKinship(x, ids, rho, recombinants = c(T,F))

    # Impossible since ids[2] is direct descendent of ids[1]:
    # p11_nr_r = twoLocusKinship(x, ids, rho, recombinants = c(F,T))
    # p11_r_r = twoLocusKinship(x, ids, rho, recombinants = c(T,T))

    k11_cs_cs = 4 * p11_nr_nr/(1-rho)^2
    if(rho > 0)
      k11_tr_cs = 4 * p11_r_nr/(rho*(1-rho))
    else k11_tr_cs = 0
    # print(k11_cs_cs)
    k11 = k11_cs_cs + k11_tr_cs
    return(k11)
  }

  stop2("Sorry, only directly linear relationships are implemented so far.")

  id1 = ids[1]
  id2 = ids[2]
  fa1 = father(x, id1)
  fa2 = father(x, id2)
  mo1 = mother(x, id1)
  mo2 = mother(x, id2)

  fou = founders(x)
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


  pp = list(c(fa1, fa2), c(fa1, mo2),c(mo1, fa2), c(mo1, mo2))
  a = matrix(NA, nrow=4, ncol=4)
  for(i in 1:4) for(j in 1:4) {
    if(i == j) {
      pi = pp[[i]]
      pj = pp[[5-i]]
      a[i,j] = phi11[pi[1], pi[2]] * phi00[pj[1], pj[2]]
    }
    else {
      pi = pp[[i]]
      pj = pp[[j]]
      a[i,j] = phi10[pi[1], pi[2]] * phi01[pj[1], pj[2]]
    }
  }

  phi11[fa1, fa2] * phi00[mo1, mo2] + phi10[fa1, fa2] * (phi01[fa1, mo2] + phi01[mo1, fa2] + phi01[mo1, mo2]) +
    phi11[fa1, mo2] * phi00[mo1, fa2] + phi10[fa1, mo2] * (phi01[fa1, fa2] + phi01[mo1, fa2] + phi01[mo1, mo2]) +
    phi11[mo1, fa2] * phi00[fa1, mo2] + phi10[mo1, fa2] * (phi01[fa1, fa2] + phi01[fa1, mo2] + phi01[mo1, mo2]) +
    phi11[mo1, mo2] * phi00[fa1, fa2] + phi10[mo1, mo2] * (phi01[fa1, fa2] + phi01[fa1, mo2] + phi01[mo1, fa2])
}
