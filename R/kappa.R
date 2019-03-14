#' IBD (kappa) coefficients
#'
#' Computes the three IBD coefficients summarising the relationship between two
#' non-inbred individuals. Both autosomal and X chromosomal versions are
#' implemented. `kappa()` is a synonym for `kappaIbd()`, but will most likely be
#' removed (since it conflicts with `base::kappa()`).
#'
#' For non-inbred individuals a and b, their autosomal IBD coefficients
#' \eqn{(\kappa0, \kappa1, \kappa2)} are defined as follows: \deqn{\kappa_i =
#' P(a and b share i alleles IBD at a random autosomal locus)}
#'
#' The autosomal kappa coefficients are computed using the method described by
#' Karigl (1981) for the condensed identity coefficients. The program first
#' checks if any of the individuals are inbred; if so, `c(NA, NA, NA)` is
#' returned. If none of the individuals are inbred, Karigl's system of equations
#' is reduced to only 3 equations: \deqn{\kappa0 + \kappa1 + \kappa2 = 1,}
#' \deqn{\kappa1 + 2*\kappa2 = 4*\phi_{ab},} \deqn{\kappa1 + 4*\kappa2 =
#' 16*\phi_{ab,ab}.} Here \eqn{\phi_{ab}} is the standard kinship coefficient,
#' and \eqn{\phi_{ab,ab}} is the generalised kinship coefficient defined by
#' Karigl (1981). The program calls [generalisedKinship22()] to compute this.
#'
#' The X chromosomal IBD coefficients are defined as in the autosomal case, with
#' the exception that \eqn{\kappa2} is undefined when at least one of the two
#' individuals is male. Hence the computation is greatly simplified when males
#' are involved. Denoting the standard kinship coefficient by \eqn{\phi}, the
#' formulas are:
#'
#' * Both male: \eqn{(\kappa0, \kappa1, \kappa2) = (1-\phi, \phi, NA)}
#'
#' * One male, one female: \eqn{(\kappa0, \kappa1, \kappa2) = (1-2*\phi, 2*\phi,
#' NA)}
#'
#' * Two females: As in the autosomal case.
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param sparse A positive integer, indicating the pedigree size limit for
#'   using sparse arrays (as implemented by the
#'   [slam](https://CRAN.R-project.org/package=slam) package) instead of
#'   ordinary arrays.
#' @param verbose A logical
#'
#' @return If `ids` has length 2: A numeric vector of length 3: \eqn{(\kappa0,
#'   \kappa1, \kappa2)}
#'
#'   If `ids` has length > 2: A data frame with one row for each pair of
#'   individuals, and 5 columns. The first two columns contain the ID labels,
#'   and columns 3-5 contain the IBD coefficients.
#'
#'   For pairs involving inbred individuals (inbred *females* in the X version)
#'   all three coefficients are reported as NA. Furthermore, the X chromosomal
#'   \eqn{kappa2} is NA whenever at least one of the two individuals is male.
#'
#'
#' @references G. Karigl (1981). _A recursive algorithm for the calculation of
#'   identity coefficients_ Annals of Human Genetics, vol. 45.
#'
#' @seealso [kinship()], [condensedIdentity()]
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

#' @rdname kappaIbd
#' @export
kappa = function(x, ids, ...) {
  warning("The function `kappa()` is renamed to `kappaIbd()` in order to avoid conflict with `base::kappa()`",
          call. = FALSE)
  kappaIbd(x, ids, ...)
}

#' @rdname kappaIbd
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
