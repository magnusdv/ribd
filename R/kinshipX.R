#' @rdname kinship
#' @export
kinshipX = function(x, ids = NULL) {
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
