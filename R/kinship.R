#' Kinship coefficients
#'
#' Compute the matrix of pairwise kinship coefficients in a pedigree. Both
#' autosomal and X-chromosomal versions are supported. The pedigree founders are
#' allowed to be inbred; see [pedtools::founderInbreeding()] for how to set this
#' up, and see Examples below.
#'
#' For two (possibly equal) members A, B of a pedigree, their autosomal (resp.
#' X-chromosomal) _kinship coefficient_ is defined as the probability that a
#' random allele from A and a random allele from B, sampled at the same
#' autosomal (resp. X-chromosomal) locus, are identical by descent relative to
#' the pedigree.
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either NULL (default), or a vector of ID labels in `x`.
#' @param simplify A logical. by default TRUE. See Value.
#' @param Xchrom A logical, indicating if the autosomal (default) or
#'   X-chromosomal kinship coefficients should be computed.
#'
#' @return A symmetric N * N matrix, where N is the number of pedigree members,
#'   or `length(ids)` if this is given, containing the pairwise kinship
#'   coefficients. If `ids` has length 2, and `simplify = TRUE`, the function
#'   returns a single number.
#'
#' @seealso [inbreeding()], [kappa()]
#'
#' @examples
#' # Kinship coefficients in a nuclear family with two children
#' x = nuclearPed(2)
#' kinship(x)
#'
#' # X chromosomal kinship coefficients in the same family
#' kinship(x, Xchrom = TRUE)
#'
#' # Autosomal kinships if the mother is 100% inbred
#' founderInbreeding(x, 2) = 1
#' kinship(x)
#'
#' # Similar for X:
#' founderInbreeding(x, 2, chromType = "X") = 1
#' kinship(x, Xchrom = TRUE)

#' @export
kinship = function(x, ids = NULL, simplify = TRUE, Xchrom = FALSE) {

  hasIds = !is.null(ids)
  singlepair = length(ids) == 2

  if(is.pedList(x)) {

    if(singlepair && simplify) {
      IDS = internalID(x, ids) # Always data frame with `id`, `comp`, `int`
      comp = IDS$comp[1]
      kin = if(comp == IDS$comp[2]) kinship(x[[comp]], ids, Xchrom = Xchrom) else 0
      return(kin)
    }

    if(!hasIds)
      ids = unlist(labels(x)) # todo: remove unlist

    if(dup <- anyDuplicated.default(ids))
      stop2("ID label is not unique: ", ids[dup])

    # Initialise big matrix
    ntot = length(ids)
    kinmat = matrix(0, nrow = ntot, ncol = ntot, dimnames = list(ids, ids))

    # Fill in component blocks
    for(comp in x) {
      idsComp = .myintersect(ids, comp$ID)
      if(length(idsComp))
        kinmat[idsComp, idsComp] = kinship(comp, ids = idsComp, simplify = FALSE, Xchrom = Xchrom)
    }

    return(kinmat)
  }

  ### Connected pedigree

  # If X, delegate to X version
  if(Xchrom)
    return(.kinshipX(x, ids = ids, simplify = simplify))

  if(!is.ped(x))
    stop2("First argument must be a `ped` object or a list of such")

  # Ensure standard order of pedigree members
  standardOrder = hasParentsBeforeChildren(x)
  if(!standardOrder) {
    origOrder = x$ID
    x = parentsBeforeChildren(x)
  }

  if(hasIds)
    IDS = internalID(x, ids) # internal index after parentsBeforeCh!

  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = which(FIDX == 0)
  NONFOU = which(FIDX > 0)
  FOU_INB = x$FOUNDER_INBREEDING[["autosomal"]] %||% rep_len(0, length(FOU)) # vector with all founders, including 0's

  # If `ids` given, restrict vectors
  if(hasIds) {
    N = max(IDS)
    FOU = FOU[FOU <= N]
    NONFOU = NONFOU[NONFOU <= N]
    length(FOU_INB) = length(FOU)
  }
  else
    N = length(FIDX) # pedsize

  # Initializing the kinship matrix.
  # Diagonal entries of founders are 0.5*(1+f)
  self_kinships = rep(0, N)
  self_kinships[FOU] = 0.5 * (1 + FOU_INB)
  kins = diag(self_kinships, nrow = N, ncol = N)

  # Vector of (maximal) generation number of each ID: dp[i] = 1 + max(dp[parents])
  # Simpler & faster than kindepth(). Requires "parentsBeforeChildren".
  dp = rep(0, N)
  for(i in NONFOU)
    dp[i] = 1 + max(dp[c(FIDX[i], MIDX[i])])

  max_dp = if(hasIds) max(dp[IDS]) else max(dp)

  # Iteratively fill the kinship matrix, one generation at a time
  for (gen in seq_len(max_dp)) {
    indx = which(dp == gen)
    Mindx = MIDX[indx]
    Findx = FIDX[indx]
    kins[indx, ] = (kins[Findx, ] + kins[Mindx, ])/2
    kins[, indx] = (kins[, Findx] + kins[, Mindx])/2
    kins[cbind(indx, indx)] = (1 + kins[cbind(Findx, Mindx)])/2
  }

  if(singlepair && simplify)
    return(kins[IDS[1], IDS[2]])

  if(hasIds) {
    kins = kins[IDS, IDS, drop = FALSE]
    dimnames(kins) = list(ids, ids)
  }
  else {
    dimnames(kins) = list(x$ID, x$ID)

    # Back to original order if needed
    if(!standardOrder)
      kins = kins[origOrder, origOrder, drop = FALSE]
  }

  kins
}


.kinshipX = function(x, ids = NULL, simplify = TRUE) {

  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  if(any(x$SEX == 0))
    stop2("Cannot compute X-kinship in pedigrees with members of unknown sex: ",
          labels(x)[x$SEX == 0])

  # Ensure standard order of pedigree members
  standardOrder = hasParentsBeforeChildren(x)
  if(!standardOrder) {
    origOrder = labels(x)
    x = parentsBeforeChildren(x)
  }

  hasIds <- !is.null(ids)
  singlepair = length(ids) == 2
  if(hasIds)
    IDS = internalID(x, ids) # internal index after parentsBeforeCh!

  FIDX = x$FIDX
  MIDX = x$MIDX
  SEX = x$SEX
  FOU = which(FIDX == 0)
  NONFOU = which(FIDX > 0)
  FOU_INB = founderInbreeding(x, chromType = "x") # vector with all founders, including 0's

  # If `ids` given, restrict vectors
  if(hasIds) {
    N = max(IDS)
    SEX = SEX[1:N]
    FOU = FOU[FOU <= N]
    NONFOU = NONFOU[NONFOU <= N]
    length(FOU_INB) = length(FOU)
  }
  else
    N = length(FIDX) # pedsize

  # Initializing the kinship matrix.
  # Diagonal entries of founders are 0.5*(1+f)
  self_kinships = rep(0, N)
  self_kinships[FOU] = ifelse(SEX[FOU] == 1, 1, 0.5 * (1 + FOU_INB))

  kins = diag(self_kinships, nrow = N, ncol = N)

  # Vector of (maximal) generation number of each ID: dp[i] = 1 + max(dp[parents])
  # Simpler & faster than kindepth(). Requires "parentsBeforeChildren".
  dp = rep(0, N)
  for(i in NONFOU)
    dp[i] = 1 + max(dp[c(FIDX[i], MIDX[i])])

  max_dp = if(singlepair) max(dp[IDS]) else max(dp)

  # Iteratively fill the kinship matrix, one generation at a time
  for (gen in seq_len(max_dp)) {

    # males
    indx_mal = which(dp == gen & SEX == 1)
    Mindx_mal = MIDX[indx_mal]
    kins[indx_mal, ] = kins[Mindx_mal, ]
    kins[, indx_mal] = kins[, Mindx_mal]
    kins[cbind(indx_mal, indx_mal)] = 1

    # females
    indx_fem = which(dp == gen & SEX == 2)
    Mindx_fem = MIDX[indx_fem]
    Findx_fem = FIDX[indx_fem]
    kins[indx_fem, ] = (kins[Findx_fem, ] + kins[Mindx_fem, ])/2
    kins[, indx_fem] = (kins[, Findx_fem] + kins[, Mindx_fem])/2
    kins[cbind(indx_fem, indx_fem)] = (1 + kins[cbind(Findx_fem, Mindx_fem)])/2
  }

  if(singlepair && simplify)
    return(kins[IDS[1], IDS[2]])



#' @rdname kinship
#' @export
kinshipX = function(x, ids = NULL) {
  message("This function is deprecated. Use `kinship(..., Xchrom = TRUE)` instead.")

  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  if(any(x$SEX == 0))
    stop2("Members with unknown gender not allowed in kinshipX algorithm: ",
          labels(x)[x$SEX == 0])

  # Ensure standard order of pedigree members
  standardOrder = hasParentsBeforeChildren(x)
  if(!standardOrder) {
    origOrder = labels(x)
    x = parentsBeforeChildren(x)
  }

  if(singlepair <- !is.null(ids)) {
    if(length(ids) != 2)
      stop2("When `ids` is not NULL, it must be a vector of length 2")
    IDS = internalID(x, ids)
  if(hasIds) {
    kins = kins[IDS, IDS, drop = FALSE]
    dimnames(kins) = list(ids, ids)
  }
  else {
    dimnames(kins) = list(x$ID, x$ID)

  FIDX = x$FIDX
  MIDX = x$MIDX
  SEX = x$SEX
  FOU = founders(x, internal = TRUE)
  NONFOU = nonfounders(x, internal = TRUE)
  N = pedsize(x)

  # Vector of X inb coeffs for all founders (including those with 0)
  FOU_INB = founderInbreeding(x, chromType = "x")

  # Initializing the kinship matrix.
  # Diagonal entries of founders are 0.5*(1+f)
  self_kinships = rep(0, N)
  self_kinships[FOU] = ifelse(SEX[FOU] == 1, 1, 0.5 * (1 + FOU_INB))

  kins = diag(self_kinships, nrow = N, ncol = N)

  # Vector of (maximal) generation number of each ID: dp[i] = 1 + max(dp[parents])
  # Simpler & faster than kindepth(). Requires "parentsBeforeChildren".
  dp = rep(0, N)
  for(i in NONFOU)
    dp[i] = 1 + max(dp[c(FIDX[i], MIDX[i])])

  max_dp = if(singlepair) max(dp[IDS]) else max(dp)

  # Iteratively fill the kinship matrix, one generation at a time
  for (gen in seq_len(max_dp)) {

    # males
    indx_mal = which(dp == gen & SEX == 1)
    Mindx_mal = MIDX[indx_mal]
    kins[indx_mal, ] = kins[Mindx_mal, ]
    kins[, indx_mal] = kins[, Mindx_mal]
    kins[cbind(indx_mal, indx_mal)] = 1

    # females
    indx_fem = which(dp == gen & SEX == 2)
    Mindx_fem = MIDX[indx_fem]
    Findx_fem = FIDX[indx_fem]
    kins[indx_fem, ] = (kins[Findx_fem, ] + kins[Mindx_fem, ])/2
    kins[, indx_fem] = (kins[, Findx_fem] + kins[, Mindx_fem])/2
    kins[cbind(indx_fem, indx_fem)] = (1 + kins[cbind(Findx_fem, Mindx_fem)])/2
    # Back to original order if needed
    if(!standardOrder)
      kins = kins[origOrder, origOrder, drop = FALSE]
  }

  if(singlepair)
    return(kins[IDS[1], IDS[2]])

  labs = labels(x)
  dimnames(kins) = list(labs, labs)

  # Back to original order if needed
  if(!standardOrder)
    kins = kins[origOrder, origOrder, drop = FALSE]

  kins
}

