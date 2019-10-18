#' Condensed identity coefficients
#'
#' Computes the 9 condensed identity coefficients of pairwise relationships in a
#' pedigree. Founders of the pedigree may be inbred; use
#' [pedtools::founderInbreeding()] to set this up.
#'
#' The implementation is a modified version of Karigl's recursive algorithm
#' (1981).
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param sparse A positive integer, indicating the pedigree size limit for
#'   using sparse arrays (as implemented by the
#'   [slam](https://CRAN.R-project.org/package=slam) package) instead of
#'   ordinary arrays.
#' @param verbose A logical
#' @param checkAnswer If TRUE, and the `identity` package is installed, the
#'   result is checked against the output of [identity::identity.coefs()].
#'   (Ignored if any of the founders are inbred.)
#'
#' @return If `ids` has length 2: A vector of length 9, containing the condensed
#'   identity coefficients.
#'
#'   If `ids` has length > 2: A data frame with one row for each pair of
#'   individuals, and 11 columns. The first two columns contain the ID labels,
#'   and columns 3-11 contain the condensed identity coefficients.
#'
#' @references G. Karigl (1981). _A recursive algorithm for the calculation of
#'   identity coefficients_ Annals of Human Genetics, vol. 45.
#'
#' @seealso [kappa()], [condensedIdentityX()], [pedtools::founderInbreeding()]
#'
#' @examples
#' # One generation of full sib mating.
#' # (This is the simplest example with all 9 coefficients nonzero.)
#' x = fullSibMating(1)
#' j1 = condensedIdentity(x, ids = 5:6)
#'
#' stopifnot(all.equal(j1, c(2, 1,4, 1, 4, 1, 7, 10, 2)/32))
#'
#' # Recalculate the coefficients when the founders are 100% inbred
#' founderInbreeding(x, 1:2) = 1
#' condensedIdentity(x, ids = 5:6)
#'
#' @importFrom utils combn
#' @export
condensedIdentity = function(x, ids, sparse = NA, verbose = FALSE, checkAnswer = verbose) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = "autosomal", verbose = verbose)

  M9 = matrix(c(
    1,1,1,1,1,1,1,1,1,
    2,2,2,2,1,1,1,1,1,
    2,2,1,1,2,2,1,1,1,
    4,0,2,0,2,0,2,1,0,
    8,0,4,0,2,0,2,1,0,
    8,0,2,0,4,0,2,1,0,
    16,0,4,0,4,0,2,1,0,
    4,4,2,2,2,2,1,1,1,
    16,0,4,0,4,0,4,1,0), byrow = T, ncol = 9)

  # If input is a pair of indivs, return the 9 coeffs as a numeric vector
  if(length(ids) == 2) {
    id1 = ids_int[1]; id2 = ids_int[2]
    RHS = c(
      1,
      2 * phi2(id1, id1, chromType = "autosomal", mem = mem),
      2 * phi2(id2, id2, chromType = "autosomal", mem = mem),
      4 * phi2(id1, id2, chromType = "autosomal", mem = mem),
      8 * phi3(id1, id1, id2, chromType = "autosomal", mem = mem),
      8 * phi3(id1, id2, id2, chromType = "autosomal", mem = mem),
      16 * phi4(id1, id1, id2, id2, chromType = "autosomal", mem = mem),
      4 * phi22(id1, id1, id2, id2, chromType = "autosomal", mem = mem),
      16 * phi22(id1, id2, id1, id2, chromType = "autosomal", mem = mem))

    j = solve(M9, RHS)

    if(verbose)
      printCounts(mem)

    if(checkAnswer)
      compare_with_identity(x, ids, j)

    return(j)
  }

  # More than 2 individuals: Do all unordered pairs; return data.frame.
  pairs = combn(ids_int, 2, simplify = F)

  RHS = vapply(pairs, function(p) {
    id1 = p[1]; id2 = p[2]
    c(1,
      2 * phi2(id1, id1, chromType = "autosomal", mem = mem),
      2 * phi2(id2, id2, chromType = "autosomal", mem = mem),
      4 * phi2(id1, id2, chromType = "autosomal", mem = mem),
      8 * phi3(id1, id1, id2, chromType = "autosomal", mem = mem),
      8 * phi3(id1, id2, id2, chromType = "autosomal", mem = mem),
      16 * phi4(id1, id1, id2, id2, chromType = "autosomal", mem = mem),
      4 * phi22(id1, id1, id2, id2, chromType = "autosomal", mem = mem),
      16 * phi22(id1, id2, id1, id2, chromType = "autosomal", mem = mem))
    }, FUN.VALUE = numeric(9))

  # Compute identity coefficients
  # Output is matrix with 9 rows
  j = solve(M9, RHS)

  # Build result data frame
  labs = labels(x)
  idcols = do.call(rbind, pairs)
  res = data.frame(id1 = labs[idcols[, 1]],
                   id2 = labs[idcols[, 2]],
                   t.default(j),
                   stringsAsFactors = F)
  names(res)[3:11] = paste0("D", 1:9)

  if(verbose)
    printCounts(mem)

  res
}

compare_with_identity = function(x, ids, j) {
  cat("Comparison with `identity` package: ")

  if(hasInbredFounders(x)) {
    message("NA (pedigree has inbred founders)")
    return()
  }

  if(!requireNamespace("identity", quietly = TRUE)) {
    message("NA (package `identity` is not installed)")
    return()
  }

  jj = idcoefs(x, ids)
  if(isTRUE(all.equal(j, jj)))
    message("OK!")
  else {
    message("*** MISMATCH! ***")
    cat("IDS:", ids, "\n")
    print(rbind(`identity:` = jj, `ribd:` = j))
  }
}

# Define double bracket extract/replace operators to accommodate the impossibility of zero-values (-1 used instead)
`[[.simple_sparse_array` = function(x, ...) {
  val = x[...]$v
  if(length(val) == 0) return(NA)
  if(val < 0) return(0)
  val
}

`[[<-.simple_sparse_array` = function(x, ..., value) {
  if(value == 0) value = -1
  x[...] <- value
  x
}
