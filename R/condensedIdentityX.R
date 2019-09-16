#' Identity coefficients on X
#'
#' Computes the X chromosomal condensed identity coefficients of a pairwise
#' relationship.
#'
#' The implementation is inspired by Karigl's recursive algorithm (1981) for the
#' autosomal case, modified to account for X-linked inheritance.
#'
#' The X chromosomal pairwise identity states depend on the sexes of the two
#' individuals. If both are female, the states are the same as in the autosomal
#' case. When males are involved, the two individuals have less than 4 alleles,
#' hence the states differ from the autosomal ones. However, to avoid drawing
#' (and learning) new pictures we re-use the autosomal states by using the
#' following simple rule: **Replace the single allele of any male, with a pair
#' of autozygous alleles**. In this way each X state corresponds to a unique
#' autosomal state.
#'
#' For simplicity the output always contain 9 coefficients, but with NA's in the
#' positions of undefined states (depending on the sex combination). The README
#' file on the github home page of ribd has a table making all of this clear.
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param verbose A Logical
#' @param sparse A positive integer, indicating the pedigree size limit for
#'   using sparse arrays (as implemented by the
#'   [slam](https://CRAN.R-project.org/package=slam) package) instead of
#'   ordinary arrays.
#'
#' @return If `ids` has length 2: A vector of length 9, containing the condensed
#'   identity coefficients. If any of the individuals are male, certain states
#'   are undefined, and the corresponding coefficients are NA. (See Details.)
#'
#'   If `ids` has length > 2: A data frame with one row for each pair of
#'   individuals, and 11 columns. The first two columns contain the ID labels,
#'   and columns 3-11 contain the condensed identity coefficients.
#'
#' @seealso [kinshipX()], [condensedIdentity()], [pedtools::founderInbreeding()]
#'
#' @examples
#' library(pedtools)
#'
#' x = fullSibMating(1)
#' x_sisters = swapSex(x, 5)
#' x_brothers = swapSex(x, 6)
#'
#' condensedIdentityX(x, ids = 5:6)
#' condensedIdentityX(x_sisters, ids = 5:6)
#' condensedIdentityX(x_brothers, ids = 5:6)
#'
#' @importFrom utils combn
#' @export
condensedIdentityX = function(x, ids, sparse = NA, verbose = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = "x", verbose = verbose)

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
      2 * phi2(id1, id1, chromType = "x", mem = mem),
      2 * phi2(id2, id2, chromType = "x", mem = mem),
      4 * phi2(id1, id2, chromType = "x", mem = mem),
      8 * phi3(id1, id1, id2, chromType = "x", mem = mem),
      8 * phi3(id1, id2, id2, chromType = "x", mem = mem),
      16 * phi4(id1, id1, id2, id2, chromType = "x", mem = mem),
      4 * phi22(id1, id1, id2, id2, chromType = "x", mem = mem),
      16 * phi22(id1, id2, id1, id2, chromType = "x", mem = mem))

    j = solve(M9, RHS)

    # Set NA at undefined states (when males are involved)
    sex = getSex(x, c(id1, id2))
    if(sex[1] == 1 && sex[2] == 1)
      j[3:9] = NA
    if(sex[1] == 1 && sex[2] == 2)
      j[5:9] = NA
    if(sex[1] == 2 && sex[2] == 1)
      j[c(3:4,7:9)] = NA

    if(verbose)
      printCounts(mem)

    return(j)
  }

  # More than 2 individuals: Do all unordered pairs; return data.frame.
  pairs = combn(ids_int, 2, simplify=F)

  RHS = vapply(pairs, function(p) {
    id1 = p[1]; id2 = p[2]
    c(1,
      2 * phi2(id1, id1, chromType = "x", mem = mem),
      2 * phi2(id2, id2, chromType = "x", mem = mem),
      4 * phi2(id1, id2, chromType = "x", mem = mem),
      8 * phi3(id1, id1, id2, chromType = "x", mem = mem),
      8 * phi3(id1, id2, id2, chromType = "x", mem = mem),
      16 * phi4(id1, id1, id2, id2, chromType = "x", mem = mem),
      4 * phi22(id1, id1, id2, id2, chromType = "x", mem = mem),
      16 * phi22(id1, id2, id1, id2, chromType = "x", mem = mem))
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

  # Set NA at undefined states (when males are involved)
  sex1 = getSex(x, res[,1])
  sex2 = getSex(x, res[,2])
  res[sex1 == 1 & sex2 == 1, 2 + 3:9] = NA
  res[sex1 == 1 & sex2 == 2, 2 + 5:9] = NA
  res[sex1 == 2 & sex2 == 1, 2 + c(3:4,7:9)] = NA

  if(verbose)
    printCounts(mem)

  res
}

