#' Kinship coefficients
#'
#' Compute the matrix of kinship coefficients of all members of a pedigree. Both
#' autosomal and X-chromosomal versions are supported. The pedigree founders are
#' allowed to be inbred; see [pedtools::founderInbreeding()] for how to set this
#' up, and see Examples below.
#'
#' For two (possibly equal) members A, B of a pedigree, their
#' autosomal (resp. X-chromosomal) _kinship coefficient_ is defined as the probability that
#' a random allele from A and a random allele from B, sampled at the same
#' autosomal (resp. X-chromosomal) locus, are identical by descent relative to
#' the pedigree.
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either NULL (default), or a vector of length 2, containing the IDs
#'   of two (possibly equal) members of `x`.
#' @param Xchrom A logical, indicating if the autosomal (default) or
#'   X-chromosomal kinship coefficients should be computed.
#'
#' @return If `ids = NULL`, a symmetric matrix containing all pairwise kinship
#'   coefficients in `x`. If `ids` has length 2, the function returns a single
#'   number.
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
kinship = function(x, ids = NULL, Xchrom = FALSE) {

  singlepair = !is.null(ids)
  if(singlepair && length(ids) != 2)
    stop2("When `ids` is not NULL, it must be a vector of length 2")

  if(is.pedList(x)) {

    if(singlepair) {  # Note: Here IDS is a data frame with cols id, comp, int
      IDS = internalID(x, ids)
      comp = IDS$comp[1]
      kin = if(comp == IDS$comp[2]) kinship(x[[comp]], ids, Xchrom = Xchrom) else 0
      return(kin)
    }

    # Initialise big matrix
    ntot = sum(pedsize(x))
    labs = unlist(labels(x))
    kinmat = matrix(0, nrow = ntot, ncol = ntot, dimnames = list(labs, labs))

    # Fill in component blocks
    for(comp in x) {
      idsComp = labels(comp)
      kinmat[idsComp, idsComp] = kinship(comp, ids = NULL, Xchrom = Xchrom)
    }

    return(kinmat)
  }

  ### Connected pedigree

  # If X, delegate to X version
  if(Xchrom)
    return(.kinshipX(x, ids = ids))

  if(!is.ped(x))
    stop2("First argument must be a `ped` object or a list of such")

  # Ensure standard order of pedigree members
  standardOrder = hasParentsBeforeChildren(x)
  if(!standardOrder) {
    origOrder = labels(x)
    x = parentsBeforeChildren(x)
  }

  if(singlepair)
    IDS = internalID(x, ids) # internal index after parentsBeforeCh!

  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = founders(x, internal = TRUE)
  NONFOU = nonfounders(x, internal = TRUE)
  N = pedsize(x)

  # Vector of inb coeffs for all founders (including those with 0)
  FOU_INB = founderInbreeding(x)

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

  max_dp = if(singlepair) max(dp[IDS]) else max(dp)

  # Iteratively fill the kinship matrix, one generation at a time
  for (gen in seq_len(max_dp)) {
    indx = which(dp == gen)
    Mindx = MIDX[indx]
    Findx = FIDX[indx]
    kins[indx, ] = (kins[Findx, ] + kins[Mindx, ])/2
    kins[, indx] = (kins[, Findx] + kins[, Mindx])/2
    kins[cbind(indx, indx)] = (1 + kins[cbind(Findx, Mindx)])/2
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


.kinshipX = function(x, ids = NULL) {

  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  if(any(x$SEX == 0))
    stop2("Members with unknown gender not allowed in X-kinship algorithm: ",
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
  }

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
  }

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

