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
#' @param simplify Simplify the output (to a numeric of length 9) if `ids` has
#'   length 2. Default: TRUE.
#' @param self A logical indicating if self-relationships (i.e., between a
#'   pedigree member and itself) should be included. FALSE by default.
#' @param verbose A logical
#'
#' @return If `ids` has length 2 and `simplify = TRUE`: A vector of length 9,
#'   containing the condensed identity coefficients.
#'
#'   Otherwise, a data frame with 11 columns and one row for each pair of
#'   individuals. The first two columns contain the ID labels, and columns 3-11
#'   contain the condensed identity coefficients.
#'
#' @references G. Karigl (1981). _A recursive algorithm for the calculation of
#'   identity coefficients._ Annals of Human Genetics, vol. 45.
#'
#' @seealso [kappa()], [identityCoefs()], [pedtools::founderInbreeding()]
#'
#' @examples
#' # One generation of full sib mating.
#' # (One of the simplest examples with all 9 coefficients nonzero.)
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
condensedIdentity = function(x, ids, sparse = NA, simplify = TRUE, self = FALSE, verbose = FALSE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")

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
    16,0,4,0,4,0,4,1,0), byrow = TRUE, ncol = 9)

  # If input is a pair of indivs, return the 9 coeffs as a numeric vector
  if(length(ids) == 2) {
    id1 = ids_int[1]; id2 = ids_int[2]
    RHS = c(
      1,
      2 * phi2(id1, id1, X = FALSE, mem = mem),
      2 * phi2(id2, id2, X = FALSE, mem = mem),
      4 * phi2(id1, id2, X = FALSE, mem = mem),
      8 * phi3(id1, id1, id2, X = FALSE, mem = mem),
      8 * phi3(id1, id2, id2, X = FALSE, mem = mem),
      16 * phi4(id1, id1, id2, id2, X = FALSE, mem = mem),
      4 * phi22(id1, id1, id2, id2, X = FALSE, mem = mem),
      16 * phi22(id1, id2, id1, id2, X = FALSE, mem = mem))

    j = solve(M9, RHS)

    if(verbose)
      printCounts(mem)

    if(simplify)
      return(j)
    else {
      res = data.frame(ids[1], ids[2], t.default(j))
      names(res) = c("id1", "id2", paste0("D", 1:9))
      return(res)
    }
  }

  # More than 2 individuals: Do all unordered pairs; return data.frame.
  pairs = .idPairs(ids_int, self = self, as = "integer")

  RHS = vapply(pairs, function(p) {
    id1 = p[1]; id2 = p[2]
    c(1,
      2 * phi2(id1, id1, X = FALSE, mem = mem),
      2 * phi2(id2, id2, X = FALSE, mem = mem),
      4 * phi2(id1, id2, X = FALSE, mem = mem),
      8 * phi3(id1, id1, id2, X = FALSE, mem = mem),
      8 * phi3(id1, id2, id2, X = FALSE, mem = mem),
      16 * phi4(id1, id1, id2, id2, X = FALSE, mem = mem),
      4 * phi22(id1, id1, id2, id2, X = FALSE, mem = mem),
      16 * phi22(id1, id2, id1, id2, X = FALSE, mem = mem))
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
                   stringsAsFactors = FALSE)
  names(res)[3:11] = paste0("D", 1:9)

  if(verbose)
    printCounts(mem)

  res
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

# Only used in legacy functions `condensedIdentity` `condensedIdentityX`
initialiseMemo = function(ped, ids, sparse = 20, chromType = "autosomal", verbose = FALSE) {

  chromType = match.arg(tolower(chromType), c("autosomal", "x"))
  Xchrom = chromType == "x"

  # Create memory storage
  mem = new.env()

  # Start timing
  st = Sys.time()

  FIDX = ped$FIDX
  MIDX = ped$MIDX
  SEX = ped$SEX

  # Logical matrix showing who has a common ancestor within the pedigree.
  ANC = hasCommonAncestor(ped)

  # Compute kinship matrix directly
  KIN = kinship(ped, Xchrom = Xchrom)
  INB = 2*diag(KIN) - 1
  if(Xchrom)
    INB[getSex(ped) == 1L] = 1 # on X, males are always 1

  maxId = max(ids)

  # When to use sparse arrays. Default: pedsize > 30
  if(is.na(sparse))
    sparse = 30

  if(is.numeric(sparse))
    sparse = maxId > sparse

  if(sparse) {
    if(verbose) message("Using sparse lookup tables")
    sparsarr0 = slam::simple_sparse_zero_array
    KIN3 = sparsarr0(dim = rep(maxId, 3), mode = "double")
    KIN4 = sparsarr0(dim = rep(maxId, 4), mode = "double")
    KIN22 = sparsarr0(dim = rep(maxId, 4), mode = "double")
  }
  else {
    KIN3 = array(NA_real_, dim = rep(maxId, 3))
    KIN4 = array(NA_real_, dim = rep(maxId, 4))
    KIN22 = array(NA_real_, dim = rep(maxId, 4))
  }

  # For quick look-up:
  FOU = founders(ped, internal = TRUE)
  isFounder = rep(FALSE, pedsize(ped))
  isFounder[FOU] = TRUE

  # Founder inbreeding
  # A vector of length pedsize(ped), with inb.coeffs at all founder idx,
  # and NA entries everywhere else. Enables quick look-up e.g. founderInb[a].
  #founderInb = rep(NA_real_, pedsize(ped))
  #founderInb[FOU] = founderInbreeding(ped, ids = founders(ped), chromType = chromType)

  for(i in FOU) {
    if(i > maxId) break # otherwise out of range!
    fi = INB[i]
    KIN3[i, i, i] = (1 + 3*fi)/4
    KIN4[i, i, i, i] = (1 + 7*fi)/8
    KIN22[i, i, i, i] = (1 + 3*fi)/4
  }

  mem$ANC = ANC
  mem$FIDX = FIDX
  mem$MIDX = MIDX
  mem$SEX = SEX

  mem$isFounder = isFounder
  #mem$founderInb = founderInb

  mem$KIN = KIN
  mem$KIN3 = KIN3
  mem$KIN4 = KIN4
  mem$KIN22 = KIN22
  mem$INB = INB

  # Counters
  mem$i2 = mem$i3 = mem$i4 = mem$i22 = 0
  mem$i2r = mem$i3r = mem$i4r = mem$i22r = 0

  # Start time
  mem$st = st

  mem
}

printCounts = function(mem) {
  tot_time = format(Sys.time()-mem$st, digits = 4)

  msg = glue::glue("
                   Function calls:
                      phi2: {mem$i2} (recursions: {mem$i2r}
                      phi3: {mem$i3} (recursions: {mem$i3r})
                      phi4: {mem$i4} (recursions: {mem$i4r})
                     phi22: {mem$i22} (recursions: {mem$i22r})
                   Lookups: {mem$ilook %||% 0}
                   Total time used: {tot_time}")
  message(msg)
}



########################################


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
#' following simple rule: **Replace any hemizygous male allele with a pair of
#' autozygous alleles**. In this way each X state corresponds to a unique
#' autosomal state.
#'
#' For simplicity the output always contains 9 coefficients, but with NA's in
#' the positions of undefined states (depending on the sex combination). The
#' README file on the GitHub home page of ribd has a table illustrating this.
#'
#' @inheritParams condensedIdentity
#'
#' @return If `ids` has length 2 and `simplify = TRUE`: A vector of length 9,
#'   containing the condensed identity coefficients. If any of the individuals
#'   are male, certain states are undefined, and the corresponding coefficients
#'   are NA. (See Details.)
#'
#'   Otherwise, a data frame with 11 columns and one row for each pair of
#'   individuals. The first two columns contain the ID labels, and columns 3-11
#'   contain the condensed identity coefficients.
#'
#' @seealso [kinship()], [identityCoefs()], [pedtools::founderInbreeding()]
#'
#' @examples
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
condensedIdentityX = function(x, ids, sparse = NA, simplify = TRUE, verbose = FALSE) {
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
    16,0,4,0,4,0,4,1,0), byrow = TRUE, ncol = 9)

  # If input is a pair of indivs, return the 9 coeffs as a numeric vector
  if(length(ids) == 2) {
    id1 = ids_int[1]; id2 = ids_int[2]
    RHS = c(
      1,
      2 * phi2(id1, id1, X = TRUE, mem = mem),
      2 * phi2(id2, id2, X = TRUE, mem = mem),
      4 * phi2(id1, id2, X = TRUE, mem = mem),
      8 * phi3(id1, id1, id2, X = TRUE, mem = mem),
      8 * phi3(id1, id2, id2, X = TRUE, mem = mem),
      16 * phi4(id1, id1, id2, id2, X = TRUE, mem = mem),
      4 * phi22(id1, id1, id2, id2, X = TRUE, mem = mem),
      16 * phi22(id1, id2, id1, id2, X = TRUE, mem = mem))

    j = solve(M9, RHS)

    # Set NA at undefined states (when males are involved)
    sex = getSex(x, ids)
    if(sex[1] == 1 && sex[2] == 1)
      j[3:9] = NA
    if(sex[1] == 1 && sex[2] == 2)
      j[5:9] = NA
    if(sex[1] == 2 && sex[2] == 1)
      j[c(3:4,7:9)] = NA

    if(verbose)
      printCounts(mem)

    if(simplify)
      return(j)
    else {
      res = data.frame(ids[1], ids[2], t.default(j))
      names(res) = c("id1", "id2", paste0("D", 1:9))
      return(res)
    }
  }

  # More than 2 individuals: Do all unordered pairs; return data.frame.
  pairs = combn(ids_int, 2, simplify = FALSE)

  RHS = vapply(pairs, function(p) {
    id1 = p[1]; id2 = p[2]
    c(1,
      2 * phi2(id1, id1, X = TRUE, mem = mem),
      2 * phi2(id2, id2, X = TRUE, mem = mem),
      4 * phi2(id1, id2, X = TRUE, mem = mem),
      8 * phi3(id1, id1, id2, X = TRUE, mem = mem),
      8 * phi3(id1, id2, id2, X = TRUE, mem = mem),
      16 * phi4(id1, id1, id2, id2, X = TRUE, mem = mem),
      4 * phi22(id1, id1, id2, id2, X = TRUE, mem = mem),
      16 * phi22(id1, id2, id1, id2, X = TRUE, mem = mem))
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
                   stringsAsFactors = FALSE)
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

