#' Kinship coefficients
#'
#' Compute the matrix of kinship coefficients of all members of a pedigree. The
#' founders may be inbred; see [pedtools::founderInbreeding()] for how to set
#' this up.
#'
#' For two (possibly identical) members $A, B$ of a pedigree, their _kinship
#' coefficient_ is defined as the probability that random alleles sampled from
#' $A$ and $B$ at the same locus, are identical by descent relative to the
#' pedigree.
#'
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#'
#' @return A symmetric matrix containing all pairwise kinship coefficients in
#'   `x`.
#'
#' @seealso [kinship()]
#' @examples
#' # Kinship coefficients in a nuclear family with two children
#' x = nuclearPed(2)
#' kinship(x)
#'
#' # Recaluclate if the father is 100% inbred
#' founderInbreeding(x, 1) = 1
#' kinship(x)
#'
#' @export
kinship = function(x) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = founders(x, internal=TRUE)
  NONFOU = nonfounders(x, internal=TRUE)

  # Vector of inb coeffs for all founders (including those with 0)
  FOU_INB = founderInbreeding(x, ids=founders(x))

  # Initializing the kinship matrix.
  # Diagonal entries of founders are 0.5*(1+f)
  self_kinships = rep(0, pedsize(x))
  self_kinships[FOU] = 0.5 * (1 + FOU_INB)

  kins = diag(self_kinships)

  # Vector of (maximal) generation number of each ID: dp[i] = 1 + max(dp[parents])
  # Gives same output as kinsip2::kindepth(), but simpler & faster.
  dp = rep(0, pedsize(x))
  for(i in NONFOU)
    dp[i] = 1 + max(dp[c(FIDX[i], MIDX[i])])

  # Iteratively fill the kinship matrix, one generation at a time
  for (gen in 1:max(dp)) {
    indx = which(dp == gen)
    Mindx = MIDX[indx]
    Findx = FIDX[indx]
    kins[indx, ] = (kins[Findx, ] + kins[Mindx, ])/2
    kins[, indx] = (kins[, Findx] + kins[, Mindx])/2
    kins[cbind(indx, indx)] = (1 + kins[cbind(Findx, Mindx)])/2
  }

  labs = labels(x)
  dimnames(kins) = list(labs, labs)
  kins
}