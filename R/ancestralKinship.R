# Implements the algorithm described in Thompson (1983), and also in
# Thompson & Morgan (1989) and the review by Thompson (1998).
# Note: The implementation below uses `ancs` as proxy for the ancestral set `S`.
# This makes everything somewhat blurry, and may work as intended only when ancs
# is a single founder, and S consists of a single allele.
# Moreover, none of the add-ons for "lethal allele" of Thompson & Morgan are
# considered here.
ancestralKinship = function(x, ids = NULL, ancs, simplify = TRUE) {

  if(!is.ped(x))
    stop2("First argument must be a `ped` object or a list of such")

  # Ensure standard order of pedigree members
  standardOrder = hasParentsBeforeChildren(x)
  if(!standardOrder) {
    origOrder = x$ID
    x = parentsBeforeChildren(x)
  }

  hasIds = !is.null(ids)
  singlepair = length(ids) == 2
  if(hasIds)
    IDS = internalID(x, ids) # internal index after parentsBeforeCh!

  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = which(FIDX == 0)
  NONFOU = which(FIDX > 0)
  N = length(FIDX) # pedsize

  if(!is.null(x$FOUNDER_INBREEDING))
    stop2("Ancestral kinship not implemented for inbred founders")

  # If `ids` given, restrict vectors
  if(hasIds) {
    N = max(IDS)
    FOU = FOU[FOU <= N]
    NONFOU = NONFOU[NONFOU <= N]
  }

  # Vector of (maximal) generation number of each ID: dp[i] = 1 + max(dp[parents])
  # Simpler & faster than kindepth(). Requires "parentsBeforeChildren".
  dp = rep(0, N)
  for(i in NONFOU)
    dp[i] = 1 + max(dp[c(FIDX[i], MIDX[i])])

  max_dp = if(hasIds) max(dp[IDS]) else max(dp)

  S = internalID(x, ancs)

  # Initialise g1 vector and g2 matrix
  g1 = g2diag = numeric(N)
  g1[FOU] = ifelse(FOU %in% S, 0.5, 0)

  g2diag[FOU] = ifelse(FOU %in% S, 0.25, 0)
  g2 = diag(g2diag, nrow = N, ncol = N)

  # Iteratively fill the kinship matrix, one generation at a time
  for (gen in seq_len(max_dp)) {
    indx = which(dp == gen)
    Mindx = MIDX[indx]
    Findx = FIDX[indx]
    g1[indx] = (g1[Mindx] + g1[Findx])/2

    g2[indx, ] = (g2[Findx, ] + g2[Mindx, ])/2
    g2[, indx] = (g2[, Findx] + g2[, Mindx])/2
    g2[cbind(indx, indx)] = (g1[indx] + g2[cbind(Findx, Mindx)])/2
  }

  if(singlepair && simplify)
    return(g2[IDS[1], IDS[2]])

  if(hasIds) {
    g2 = g2[IDS, IDS, drop = FALSE]
    dimnames(g2) = list(ids, ids)
  }
  else {
    dimnames(g2) = list(x$ID, x$ID)

    # Back to original order if needed
    if(!standardOrder)
      g2 = g2[origOrder, origOrder, drop = FALSE]
  }

  g2
}

