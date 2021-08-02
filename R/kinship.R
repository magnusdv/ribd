#' Kinship coefficients
#'
#' Compute the matrix of kinship coefficients (autosomal or X) of all members of
#' a pedigree. The founders may be inbred; see [pedtools::founderInbreeding()]
#' for how to set this up.
#'
#' For two (not necessarily distinct) members A, B of a pedigree, their
#' autosomal (resp. X) _kinship coefficient_ is defined as the probability that
#' random alleles sampled from A and B at the same autosomal (resp. X) locus,
#' are identical by descent relative to the pedigree.
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either a character of length 2, or NULL. In the former case, it
#'   must contain the ID labels of two members of `x`, and the function will
#'   return their kinship coefficient as a single number. If `ids` is NULL (this
#'   is the default), the output is the complete kinship matrix.
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
#' kinshipX(x)
#'
#' # Recalculate the autosomal kinships if the father is 100% inbred
#' founderInbreeding(x, 1) = 1
#' kinship(x)
#'
#' @export
kinship = function(x, ids = NULL) {
  if(is.pedList(x)) {

    if(length(ids) == 2) {
      compNr = getComponent(x, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
      if(compNr[1] != compNr[2])
        return(0)
      else
        return(kinship(x[[compNr[1]]], ids))
    }

    # Initialise big matrix
    ntot = sum(pedsize(x))
    labs = unlist(labels(x))
    kinmat = matrix(0, nrow = ntot, ncol = ntot, dimnames = list(labs, labs))

    # Fill in component blocks
    for(comp in x) {
      idsComp = labels(comp)
      kinmat[idsComp, idsComp] = kinship(comp, ids = NULL)
    }

    return(kinmat)
  }
  else if(!is.ped(x))
    stop2("First argument must be a `ped` object or a list of such")

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

  if(singlepair) return(kins[IDS[1], IDS[2]])

  labs = labels(x)
  dimnames(kins) = list(labs, labs)

  # Back to original order if needed
  if(!standardOrder)
    kins = kins[origOrder, origOrder, drop = FALSE]

  kins
}
