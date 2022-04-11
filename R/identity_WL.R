#####################################################################
# Computation of condensed identity coefficients
# Algorithm by Lange & Sinsheimer (1992), based on generalised kinship coefficients by Weeks & Lange (1988)
#
# Note: For detailed coefficients, see `identity_LS()`
#####################################################################

identity_WL = function(x, ids, Xchrom = FALSE, self = FALSE, mem = NULL, verbose = FALSE) {

  x = prepPed(x)
  pairs = prepIds(x, ids, self = self)

  if(is.null(mem))
    mem = memoIdentity(x, method = "WL")

  # Recursion wrapper to simplify typing
  recu = function(..., debug = FALSE) {
    recurse_WL(kinpat(x, list(...)), X = Xchrom, mem = mem, debug = debug)
  }

  # Function for computing all Phi's for specific pair
  Phi = function(id1, id2) {
    c(S1 = recu(c(id1, id1, id2, id2)),
      S2 = recu(c(id1, id1), c(id2, id2)),
      S3 = recu(c(id1, id1, id2), id2),
      S4 = recu(c(id1, id1), id2, id2),
      S5 = recu(c(id1, id2, id2), id1),
      S6 = recu(id1, id1, c(id2, id2)),
      S7 = recu(c(id1, id2), c(id1, id2)),
      S8 = recu(c(id1, id2), id1, id2),
      S9 = recu(id1, id1, id2, id2))
  }

  # Compute Phi for each pair
  PhiMat = vapply(pairs, function(p) Phi(p[1], p[2]), FUN.VALUE = numeric(9))

  # Phi -> Psi (see Lange, chapter 5.5)
  PsiMat = PhiMat * c(1, 1, 2, 1, 2, 1, 2, 4, 1)

  # Psi -> Delta triangular matrix
  B = matrix(nrow = 9, ncol = 9, data = c(
    1, 0, -.5, 0,  -.5, 0,  .5, .25, 0,
    0, 1, -.5, -1, -.5, -1, .5, .75, 1,
    0, 0, 2, 0, 0, 0, -2, -1, 0,
    0, 0, 0, 2, 0, 0, 0, -1, -2,
    0, 0, 0, 0, 2, 0, -2, -1, 0,
    0, 0, 0, 0, 0, 2, 0, -1, -2,
    0, 0, 0, 0, 0, 0, 4, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 4, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 4))

  # Delta
  j = t.default(PsiMat) %*% B

  if(verbose)
    printMemInfo(mem)

  # Output data frame
  jmat2df(j, pairs)
}

gKinship_WL = function(x, kp, Xchrom = FALSE, mem, debug = FALSE) {
  if(Xchrom)
    stop2("X is not implemented in method 'WL' yet")

  recurse_WL(kp, X = Xchrom, mem = mem, debug = debug)
}


# Recursion method of Weeks & Lange
recurse_WL = function(kp, X = FALSE, mem = NULL, debug = FALSE, indent = 0) {
  if(debug)
    cat(strrep(" ", indent), kinpat2string(kp), "\n", sep = "")
  mem$i = mem$i + 1

  kp = kinpatReduce(kp)
  L = lengths(kp)

  # B0: Trivial pattern
  if(sum(L) <= 1)  {
    mem$B0 = mem$B0 + 1
    return(debugReturn(1, debug = debug, indent = indent, comment = " (B0)"))
  }

  # B2: Any group with 2 unrelated indivs?
  uniqList = lapply(kp, unique.default)
  if(boundaryB2(kp, mem$REL, uniqList)) {
    mem$B2 = mem$B2 + 1
    return(debugReturn(0, debug = debug, indent = indent, comment = " (B2)"))
  }

  # B1: Anyone in >2 groups
  uniqVec = unlist(uniqList, use.names = FALSE)
  tab = tabulate(uniqVec)
  if(any(tab > 2)) {
    mem$B1 = mem$B1 + 1
    return(debugReturn(0, debug = debug, indent = indent, comment = " (B1)"))
  }

  # Boundary 3: All founders. (Extended to account for founder inbreeding.)
  if(all(mem$isFounder[uniqVec])) {
    mem$B3 = mem$B3 + 1
    INB = mem$INB
    tabtot = tabulate(unlist(kp, use.names = FALSE))
    res = 1
    if(length(one <- which(tab == 1)))  # those in 1 group: may be autozygous
      res = res * prod(INB[one] * 1 + (1 - INB[one]) * .5^(tabtot[one] - 1))
    if(length(two <- which(tab == 2)))  # those in 2 groups: must be non-autozygous
      res = res * prod((1 - INB[two]) * .5^(tabtot[two] - 1))
    # Original: res = 1/2^(sum(L) - length(unique.default(uniqVec)))
    return(debugReturn(res, debug = debug, indent = indent, comment = " (B3)"))
  }

  # Wrapper of recurse_WL to save typing in the recursions
  recu = function(k) recurse_WL(k, X = X, mem = mem, debug = debug, indent = indent +2)

  # If 2 (for simplicity) unrelated groups: Factorise
  if(length(kp) == 2) {
    u1 = uniqList[[1]]; u2 = uniqList[[2]]
    interpairs = cbind(rep(u1, each = length(u2)), rep(u2, length(u1)))
    if(!any(mem$REL[interpairs])) { # any pair not related?
      res = recu(`[[<-`(kp, 2, NULL)) * recu(`[[<-`(kp, 1, NULL))
      return(debugReturn(res, debug = debug, indent = indent))
    }
  }

  # Sort
  kp = kinpatSort(kp)

  # Lookup in array; compute if necessary.
  kinStr = paste(kp, collapse = ", ")
  val = mem$MEM[[kinStr]]

  if(!is.null(val)) {
    mem$ilook = mem$ilook + 1
    return(debugReturn(val, debug = debug, indent = indent, comment = " (lookup)"))
  }

  # Prepare recursion
  mem$irec = mem$irec + 1

  pivot = kp[[1]][1]
  fa = mem$FIDX[pivot]
  mo = mem$MIDX[pivot]

  # Number of times the pivot occurs in first block (s) and second block (t)
  s = sum(kp[[1]] == pivot)
  t = if(length(kp) > 1) sum(kp[[2]] == pivot) else 0

  # Recurrence rules 1 (s = 1) and 2
  if(t == 0) {
    A1 = .5^s
    res =
      A1 * recu(kinpatReplace(kp, id = pivot, rep1 = fa)) +
      A1 * recu(kinpatReplace(kp, id = pivot, rep1 = mo))

    if(s > 1) {
      B1 = (1 - 2 * (.5)^s)
      res = res +
        B1 * recu(kinpatReplace(kp, id = pivot, rep1 = c(fa, mo)))
    }
  }

  # Recurrence rule 3
  if(t > 0) {
    A2 = .5^(s + t)
    res =
      A2 * recu(kinpatReplace(kp, id = pivot, rep1 = fa, rep2 = mo)) +
      A2 * recu(kinpatReplace(kp, id = pivot, rep1 = mo, rep2 = fa))
  }

  mem$MEM[[kinStr]] = res
  return(debugReturn(res, debug = debug, indent = indent))
}


# W&L boundary 2: Group with 2 unrelated --> 0
boundaryB2 = function(kp, REL, uniqList = lapply(kp, unique.default)) {
  for(s in uniqList) {
    if(length(s) < 2)
      next

    pairMat = comb2(s, vec = TRUE)
    if(any(!REL[pairMat]))
      return(TRUE)
  }
  FALSE
}


