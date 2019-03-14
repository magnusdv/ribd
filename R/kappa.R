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

#' @export
kappaIbd = function(x, ids, sparse = NA, verbose = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = "autosomal", verbose = verbose)
  KIN2 = mem$KIN2
  INB = 2*diag(KIN2) - 1 # inbreeding coeffs

  if(verbose && any(INB[ids] > .Machine$double.eps)) {
    message("Warning: kappas involving inbred individuals are undefined:")
    inbred = ids[INB[ids] > .Machine$double.eps]
    for(id in inbred)
      message(sprintf("  %s: f = %f", id, INB[id]))
    cat("\n")
  }

  # All unordered pairs
  pairs = combn(ids_int, 2, simplify=F)

  # System of equations:
  # k0 + k1 +   k2 = 1
  #      k1 + 2*k2 = 4*phi_{ab}
  #      k1 + 4*k2 = 16*phi_{ab,ab}
  # ==> Solve for k0, k1, k2

  # Compute kappa coefficients
  kappas = vapply(pairs, function(p) {
    if(any(INB[p] > .Machine$double.eps)) {
      c(NA_real_, NA_real_, NA_real_)
    }
    else {
      id1 = p[1]; id2 = p[2]
      u = KIN2[[id1, id2]]
      v = phi22(id1, id2, id1, id2, chromType = "autosomal", mem = mem)
      c(1 - 6*u + 8*v, 8*u - 16*v, 8*v - 2*u)
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

#' @export
kappaIbdX = function(x, ids, sparse = NA, verbose = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

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
    cat("\n")
  }

  # All unordered pairs
  pairs = combn(ids_int, 2, simplify=F)

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

    id1 = p[1]; id2 = p[2]
    u = KIN2[[id1, id2]]
    switch(sum(SEX[p] == 2) + 1,
           c(1 - u, u, NA), # both males
           c(1 - 2*u, 2*u, NA), # one male, one female
           { # both female
             v = phi22(id1, id2, id1, id2, chromType = "x", mem = mem)
             c(1 - 6*u + 8*v, 8*u - 16*v, 8*v - 2*u)
           })
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
