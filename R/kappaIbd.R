#'IBD (kappa) coefficients
#'
#'Computes the three IBD coefficients summarising the relationship between two
#'non-inbred individuals. Both autosomal and X chromosomal versions are
#'implemented.
#'
#'For non-inbred individuals a and b, their autosomal IBD coefficients
#'\eqn{(\kappa0, \kappa1, \kappa2)} are defined as follows: \deqn{\kappa_i = P(a
#'and b share i alleles IBD at a random autosomal locus)}
#'
#'The autosomal kappa coefficients are computed from the kinship coefficients.
#'When a and b are both nonfounders, the following formulas are well-known:
#'
#'* \eqn{\kappa2 = \phi_MM * \phi_FF + \phi_MF * \phi_FM}
#'
#'* \eqn{\kappa1 = 4 * \phi_ab - 2 * \kappa2}
#'
#'* \eqn{\kappa0 = 1 - \kappa1 - \kappa2}
#'
#'Here \eqn{\phi_MM} denotes the kinship coefficient between the mothers of a
#'and b, and so on. If either a or b is a founder, then \eqn{\kappa2 = 0}, while
#'the other two formulas remain as above.
#'
#'The X chromosomal IBD coefficients are defined as in the autosomal case, with
#'the exception that \eqn{\kappa2} is undefined when at least one of the two
#'individuals is male. Hence the computation is greatly simplified when males
#'are involved. Denoting the standard kinship coefficient by \eqn{\phi}, the
#'formulas are:
#'
#'* Both male: \eqn{(\kappa0, \kappa1, \kappa2) = (1-\phi, \phi, NA)}
#'
#'* One male, one female: \eqn{(\kappa0, \kappa1, \kappa2) = (1-2*\phi, 2*\phi,
#'NA)}
#'
#'* Two females: As in the autosomal case.
#'
#'@param x A pedigree in the form of a [`pedtools::ped`] object.
#'@param ids A character (or coercible to character) containing ID labels of two
#'  or more pedigree members.
#'@param inbredAction An integer telling the program what to do if either of the
#'  `ids` individuals are inbred. Possible values are: 0 = do nothing; 1 = print
#'  a warning message (default); 2 = raise an error. In the first two cases
#'  the coefficients are reported as `NA`.
#'@param sparse A positive integer, indicating the pedigree size limit for using
#'  sparse arrays (as implemented by the
#'  [slam](https://CRAN.R-project.org/package=slam) package) instead of ordinary
#'  arrays.
#'@param verbose A logical.
#'@param ... Further arguments.
#'
#'@return If `ids` has length 2: A numeric vector of length 3: \eqn{(\kappa0,
#'  \kappa1, \kappa2)}.
#'
#'  If `ids` has length > 2: A data frame with one row for each pair of
#'  individuals, and 5 columns. The first two columns contain the ID labels, and
#'  columns 3-5 contain the IBD coefficients.
#'
#'  Unless `inbredAction = 2`, the coefficients of pairs involving inbred
#'  individuals (inbred *females* in the X version) are reported as NA.
#'  Furthermore, the X chromosomal \eqn{\kappa2} is NA whenever at least one of
#'  the two individuals is male.
#'
#' @seealso [kinship()], [condensedIdentity()]
#' @examples
#' # Siblings
#' x = nuclearPed(2)
#' k = kappaIBD(x, 3:4)
#' stopifnot(identical(k, c(.25, .5, .25)))
#'
#' # Quad half first cousins
#' x = quadHalfFirstCousins()
#' k = kappaIBD(x, leaves(x))
#' stopifnot(identical(k, c(17/32, 14/32, 1/32)))
#'
#' # Paternal half brothers with 100% inbred father
#' # Genetically indistinguishable from an (outbred) father-son relationship
#' x = halfSibPed()
#' founderInbreeding(x, 1) = 1
#'
#' k = kappaIBD(x, 4:5)
#' stopifnot(identical(k, c(0, 1, 0)))
#'
#'@export
kappaIBD = function(x, ids = labels(x), inbredAction = 1) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  labs = labels(x)
  ids = as.character(ids)
  if(!all(ids %in% labs))
    stop2("Unknown ID label: ", setdiff(ids, labs))
  if(length(ids) < 2)
    stop2("At least two ID labels must be indicated")
  if(dup <- anyDuplicated.default(ids))
    stop2("Duplicated ID label: ", ids[dup])

  KIN = kinship(x)
  INB = 2*diag(KIN) - 1 # inbreeding coeffs

  isInbred = INB[ids] > .Machine$double.eps
  if(any(isInbred) && inbredAction > 0) {
    msg = paste0(c(sprintf(" Individual '%s' is inbred (f = %g)", ids[isInbred], INB[ids[isInbred]]),
                  "Kappa coefficients are only defined for non-inbred individuals."), collapse = "\n")
    switch(inbredAction, {message("Warning:"); message(msg)}, stop2(msg))
  }

  # Build result data frame
  pairs = t.default(combn(ids, 2))
  res = data.frame(id1 = pairs[, 1], id2 = pairs[, 2],
                   kappa0 = NA_real_, kappa1 = NA_real_, kappa2 = NA_real_,
                   stringsAsFactors = F)

  # Rows that needs computing
  founder_rows = res$id1 %in% founders(x) | res$id2 %in% founders(x)
  noninbred_rows = INB[res$id1] < .Machine$double.eps &
    INB[res$id2] < .Machine$double.eps

  # Noninbred pairs involving founder(s)
  if(any(noninb_fou_rows <- noninbred_rows & founder_rows)) {
    k1_fou = 4*KIN[pairs[noninb_fou_rows, , drop = F]]
    res[noninb_fou_rows, 3:5] = cbind(1 - k1_fou, k1_fou, 0)
  }

  # Noninbred nonfounders
  if(any(nn_rows <- noninbred_rows & !founder_rows)) {
    id1.nn = res$id1[nn_rows]
    id2.nn = res$id2[nn_rows]
    F1 = father(x, id1.nn)
    M1 = mother(x, id1.nn)
    F2 = father(x, id2.nn)
    M2 = mother(x, id2.nn)

    k2 = KIN[cbind(F1, F2)]*KIN[cbind(M1, M2)] + KIN[cbind(F1, M2)]*KIN[cbind(M1, F2)]
    k1 = 4*KIN[pairs[nn_rows, , drop = F]] - 2*k2
    k0 = 1 - k1 - k2
    res[nn_rows, 3:5] = cbind(k0, k1, k2)
  }

  if(length(ids) == 2)
    return(as.numeric(res[1, 3:5]))

  res
}


#' @rdname kappaIBD
#' @export
kappaIbdX = function(x, ids, sparse = NA, verbose = FALSE) {
  # TODO!!: Simplify this as in the autosomal case
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = "x", verbose = verbose)
  KIN2 = mem$KIN2
  INB = 2*diag(KIN2) - 1 # inbreeding coeffs
  SEX = mem$SEX

  if(verbose && any(INB[SEX == 2] > .Machine$double.eps)) {
    message("Warning: X-kappas involving inbred females are undefined:")
    inbred = females(x)[INB[SEX == 2] > .Machine$double.eps]
    for(id in inbred)
      message(sprintf("  %s: f = %f", id, INB[id]))
    message("")
  }

  # All unordered pairs
  pairs = combn(ids_int, 2, simplify = F)

  # System of equations:
  # k0 + k1 +   k2 = 1
  #      k1 + 2*k2 = 4*phi_{ab}
  #      k1 + 4*k2 = 16*phi_{ab,ab}
  # ==> Solve for k0, k1, k2

  # Compute kappa coefficients
  kappas = vapply(pairs, function(p) {

    # If the pair includes an inbred female, return NA's
    if(any(SEX[p] == 2 & INB[p] > .Machine$double.eps))
      c(NA_real_, NA_real_, NA_real_)
    else {
      id1 = p[1]; id2 = p[2]
      u = KIN2[[id1, id2]]
      switch(sum(SEX[p] == 2) + 1,
             c(1 - u, u, NA), # both males
             c(1 - 2*u, 2*u, NA), # one male, one female
             { # both female
               v = phi22(id1, id2, id1, id2, chromType = "x", mem = mem)
               c(1 - 6*u + 8*v, 8*u - 16*v, 8*v - 2*u)
             })
    }
  }, FUN.VALUE = numeric(3))

  if(verbose)
    printCounts(mem)

  if(length(ids) == 2)
    return(kappas[,1])

  # Build result data frame
  labs = labels(x)
  idcols = do.call(rbind, pairs)
  res = data.frame(id1 = labs[idcols[, 1]],
                   id2 = labs[idcols[, 2]],
                   t.default(kappas),
                   stringsAsFactors = F)
  names(res)[3:5] = paste0("kappa", 0:2)

  res
}
