#' Jacquard coefficients
#'
#' This function calculates the pairwise identity coefficients defined by
#' Jacquard (1974). Unlike [condensedIdentity()] (and
#' `identity::identity.coefs()`), this function also computes the 15 *detailed*
#' identity coefficients. If `detailed = FALSE` the output is equivalent to
#' [condensedIdentity()], which is usually faster.
#'
#' Karigl (1981) gave the first recursive algorithm for the 9 condensed Jacqaurd
#' coefficients, based on certain *generalised kinship coefficients*. Karigl's
#' method is implemented in [condensedIdentity()].
#'
#' Weeks & Lange (1988) suggested a broader and more natural generalisation of
#' kinship coefficients, and a recursive algorithm for computing them. This is
#' implemented in [generalisedKinship()].
#'
#' Lange & Sinsheimer (1992) described an even further generalisation of kinship
#' coefficients, allowing a mix of deterministic and random sampling of alleles.
#' This is also implemented in [generalisedKinship()].
#'
#' In the same paper, Lange & Sinsheimer (1992) used the new generalisations to
#' give (i) an alternative algorithm for the 9 condensed identity coefficients,
#' and (ii) an algorithm for the 15 detailed coefficients. Both of these are
#' implemented in `jacquard()`.
#'
#' Both the condensed and detailed coefficients are given in the orders used by
#' Jacquard (1974). The function `detailed2condensed` converts from detailed
#' coefficients (d1, ... d15) to condensed ones (D1, ..., D9) using the
#' following relations:
#'
#' * D1 = d1
#'
#' * D2 = d6
#'
#' * D3 = d2 + d3
#'
#' * D4 = d7
#'
#' * D5 = d4 + d5
#'
#' * D6 = d8
#'
#' * D7 = d9 + d12
#'
#' * D8 = d10 + d11 + d13 + d14
#'
#' * D9 = d15
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object.
#' @param ids A vector of two ID labels.
#' @param detailed A logical. If FALSE (default), the 9 condensed coefficients
#'   are computed; otherwise the 15 detailed identity coefficients.
#' @param named A logical, indicating if the output vector should be named.
#' @param verbose A logical.
#' @param d A numerical vector of length 15.
#'
#' @return A numeric vector of length 9 or 15 (depending on `detailed`). The
#'   coefficients are given in the order used by Jacquard (1974).
#'
#' @seealso [condensedIdentity()], [idcoefs()]
#' @references
#'
#' * Jacquard, A. (1974). The Genetic Structure of Populations. Springer.
#'
#' * Karigl, G. (1981). A recursive algorithm for the calculation of identity
#' coefficients. Ann. Hum. Genet.
#'
#' * Weeks, D.E. & Lange, K. (1988). The affected-pedigree-member method of
#' linkage analysis. Am. J. Hum. Genet
#'
#' * Lange, K. & Sinsheimer, J.s. (1992). Calculation of genetic identity
#' coefficients. Ann. Hum. Genet.
#'
#' @examples
#' x = fullSibMating(1)
#' ids = 5:6
#'
#' # Condensed coefficient
#' j = jacquard(x, ids)
#' stopifnot(identical(j, condensedIdentity(x, ids)))
#'
#' # Detailed coefficients
#' jdet = jacquard(x, ids, detailed = TRUE)
#' jdet
#' stopifnot(identical(j, detailed2condensed(jdet)))
#'
#' @export
jacquard = function(x, ids, detailed = FALSE, named = FALSE, verbose = FALSE) {

  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  if(detailed) {
    for(id in intersect(ids, founders(x)))
      x = addParents(x, id, verbose = FALSE)
  }

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  x = foundersFirst(x)

  # Setup memoisation
  mem = initialiseGKMemo(x, counters = c("i", "itriv", "iimp", "ifound", "ilook", "irec"))

  id1 = ids[1]
  id2 = ids[2]

  if(!detailed) {
    states = list(
      S1 = kinPattern(x, list(c(id1, id1, id2, id2))),
      S2 = kinPattern(x, list(c(id1, id1), c(id2, id2))),
      S3 = kinPattern(x, list(c(id1, id1, id2), id2)),
      S4 = kinPattern(x, list(c(id1, id1), id2, id2)),
      S5 = kinPattern(x, list(c(id1, id2, id2), id1)),
      S6 = kinPattern(x, list(id1, id1, c(id2, id2))),
      S7 = kinPattern(x, list(c(id1, id2), c(id1, id2))),
      S8 = kinPattern(x, list(c(id1, id2), id1, id2)),
      S9 = kinPattern(x, list(id1, id1, id2, id2))
    )

    Phi = vapply(states, function(s) genKin(s, mem = mem, indent = NA), FUN.VALUE = 0)

    # See Lange, chapter 5.5
    Psi = Phi * c(1, 1, 2, 1, 2, 1, 2, 4, 1)

    coeffs = c(
      D1 = Psi %*% c(1, 0, -.5, 0,  -.5, 0,  .5, .25, 0),
      D2 = Psi %*% c(0, 1, -.5, -1, -.5, -1, .5, .75, 1),
      D3 = Psi %*% c(0, 0, 2, 0, 0, 0, -2, -1, 0),
      D4 = Psi %*% c(0, 0, 0, 2, 0, 0, 0, -1, -2),
      D5 = Psi %*% c(0, 0, 0, 0, 2, 0, -2, -1, 0),
      D6 = Psi %*% c(0, 0, 0, 0, 0, 2, 0, -1, -2),
      D7 = 4 * Psi[[7]],
      D8 = 4 * Psi[[8]],
      D9 = 4 * Psi[[9]]
    )

    names(coeffs) = if(named) paste0("D", 1:9) else NULL
  }
  else {
    # Detailed identity states (following Figure 6.2 (page 105), Jacquard 1974)
    states = list(
      s1 = kinPattern(x, list(c(p=id1, m=id1, p=id2, m=id2))),
      s2 = kinPattern(x, list(c(p=id1, m=id1, p=id2), c(m=id2))),
      s3 = kinPattern(x, list(c(p=id1, m=id1, m=id2), c(p=id2))),
      s4 = kinPattern(x, list(c(p=id1, p=id2, m=id2), c(m=id1))),
      s5 = kinPattern(x, list(c(m=id1, p=id2, m=id2), c(p=id1))),
      s6 = kinPattern(x, list(c(p=id1, m=id1), c(p=id2, m=id2))),
      s7 = kinPattern(x, list(c(p=id1, m=id1), c(p=id2), c(m=id2))),
      s8 = kinPattern(x, list(c(p=id1), c(m=id1), c(p=id2, m=id2))),
      s9 = kinPattern(x, list(c(p=id1, p=id2), c(m=id1, m=id2))),
      s10 = kinPattern(x, list(c(p=id1, p=id2), c(m=id1), c(m=id2))),
      s11 = kinPattern(x, list(c(p=id1), c(p=id2), c(m=id1, m=id2))),
      s12 = kinPattern(x, list(c(p=id1, m=id2), c(m=id1, p=id2))),
      s13 = kinPattern(x, list(c(p=id1, m=id2), c(m=id1), c(p=id2))),
      s14 = kinPattern(x, list(c(p=id1), c(m=id2), c(m=id1, p=id2))),
      s15 = kinPattern(x, list(c(p=id1), c(m=id1), c(m=id2), c(p=id2)))
    )

    coeffs = vapply(states, function(s) genKinDetailed(s, mem = mem, indent = NA), FUN.VALUE = 0)
    names(coeffs) = if(named) paste0("d", 1:15) else NULL
  }

  if(verbose)
    printCounts2(mem)

  coeffs
}


#' @rdname jacquard
#' @export
detailed2condensed = function(d, named = FALSE) {
  cond = c(
    d[[1]],
    d[[6]],
    d[[2]] + d[[3]],
    d[[7]],
    d[[4]] + d[[5]],
    d[[8]],
    d[[9]] + d[[12]],
    d[[10]] + d[[11]] + d[[13]] + d[[14]],
    d[[15]]
  )
  names(cond) = if(named) paste0("D", 1:9) else NULL
  cond
}
