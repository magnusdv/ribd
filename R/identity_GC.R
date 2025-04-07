#############################################
# Recursive algorithm by Garcia-Cortes (2015)
#
# Note: Condensed coeffs are calculated from
# the detailed ones
#############################################

identity_GC = function(x, ids, Xchrom = FALSE, detailed = FALSE, self = FALSE, mem = NULL, verbose = FALSE) {

  x = prepPed(x)
  pairs = prepIds(x, ids, self = self)

  if(is.null(mem))
    mem = memoIdentity(x, Xchrom = Xchrom, method = "GC")

  mem$ilook = mem$irec = mem$i = integer(4)

  # Recursion wrapper to simplify typing
  recu = function(..., debug = FALSE) {
    recurse_GC(gip(x, list(...), distinct = TRUE), X = Xchrom, mem = mem, debug = debug)
  }

  # Garcia-Cortes, equation 4
  Psi = function(id1, id2) {
    c(s1 = recu(c(p=id1, m=id1, p=id2, m=id2)),
      s2 = recu(c(p=id1, m=id1, p=id2)),
      s3 = recu(c(p=id1, m=id1, m=id2)),
      s4 = recu(c(p=id1, p=id2, m=id2)),
      s5 = recu(c(m=id1, p=id2, m=id2)),
      s6 = recu(c(p=id1, m=id1), c(p=id2, m=id2)),
      s7 = recu(c(p=id1, m=id1)),
      s8 = recu(c(p=id2, m=id2)),
      s9 = recu(c(p=id1, p=id2), c(m=id1, m=id2)),
      s10 = recu(c(p=id1, p=id2)),
      s11 = recu(c(m=id1, m=id2)),
      s12 = recu(c(p=id1, m=id2), c(m=id1, p=id2)),
      s13 = recu(c(p=id1, m=id2)),
      s14 = recu(c(m=id1, p=id2)),
      s15 = 1)
  }

  PsiMat = vapply(pairs, function(p) Psi(p[1], p[2]), numeric(15))

  # M15 = matrix(byrow = TRUE, ncol = 15, c(
  #   1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  #   1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
  #   1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
  #   1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
  #   1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
  #   1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
  #   1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,
  #   1,0,0,1,1,1,0,1,0,0,0,0,0,0,0,
  #   1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
  #   1,1,0,1,0,0,0,0,1,1,0,0,0,0,0,
  #   1,0,1,0,1,0,0,0,1,0,1,0,0,0,0,
  #   1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
  #   1,0,1,1,0,0,0,0,0,0,0,1,1,0,0,
  #   1,1,0,0,1,0,0,0,0,0,0,1,0,1,0,
  #   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))

  Minv = matrix(byrow = TRUE, nrow = 15, c(
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # d1
    -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # d2
    -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # d3
    -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # d4
    -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # d5
    -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, # d6
    2,-1,-1, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, # d7
    2, 0, 0,-1,-1,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, # d8
    -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, # d9
    2,-1, 0,-1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, # d10
    2, 0,-1, 0,-1, 0, 0, 0,-1, 0, 1, 0, 0, 0, 0, # d11
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, # d12
    2, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, # d13
    2,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0,-1, 0, 1, 0, # d14
    -6, 2, 2, 2, 2, 1,-1,-1, 1,-1,-1, 1,-1,-1, 1  # d15
  ))

  j = Minv %*% PsiMat

  if(verbose)
    printMemInfo(mem)

  # Data frame
  res = jmat2df(t.default(j), pairs)

  if(!detailed)
    res = detailed2condensed(res)

  if(Xchrom)
    res = Xmask(x, res)

  res
}


gKinship_GC = function(x, gp, Xchrom = FALSE, mem, debug = FALSE) {
  if(!isDeterministic(gp))
    stop2("Method 'GC' supports deterministic kinship patterns only")

  if(length(gp) > 1 && isDistinct(gp))
    stop2("The `GC` method requires distinct blocks. Received: ", gip2string(gp))

  recurse_GC(gp, X = Xchrom, mem = mem, debug = debug)
}

recurse_GC = function(gp, X = FALSE, mem = NULL, debug = FALSE, indent = 0) {
  if(length(gp) == 1)
    psi1G(gp, X = X, mem = mem, debug = debug)
  else if(identical(lengths(gp), c(2L, 2L)))
    psi22(gp, X = X, mem = mem, debug = debug)
  else
    stop2("Method 'GC' currently only supports patterns of either 1 group or 2 groups with two elements")
}

# Main recursion function for generalised gamete probabilities with 1 group
psi1G = function(gp, X = FALSE, mem = NULL, debug = FALSE, indent = 0) {
  if(debug)
    cat(strrep(" ", indent), gip2string(gp), "\n", sep = "")

  g = gp[[1]]
  L = length(g)

  # Quick return if trivial case
  if(L == 1)
    return(debugReturn(1, debug = debug, indent = indent, comment = " (trivial)"))


  mem$i[L] = mem$i[L] + 1L
  ids = g %/% 10

  # X males: Convert any paternal alleles to maternal. (Instead of raising error or returning NA.)
  # This simplifies things e.g. in identityCoefs, allowing autosomal code to be reused.
  if(X) {
    malePat = mem$SEX[ids] == 1 & g %% 10 == 1
    if(any(malePat))
      g[malePat] = g[malePat] + 1L
  }

  # Quick return if two alleles from same person
  if(L == 2) {
    if(g[1] == g[2])
      return(debugReturn(1, debug = debug, indent = indent, comment = " (B0)"))

    if(abs(g[1] - g[2]) == 1) #  p & m from same person
      return(debugReturn(mem$INB[[ids[1]]], debug = debug, indent = indent, comment = " (B1)"))
  }

  # Any duplicates -> remove and recurse
  if(dup <- anyDuplicated.default(g)) {
    gp[[1]] = g[-dup]
    return(psi1G(gp, X = X, mem = mem, debug = debug, indent = indent + 2))
  }

  # Any pair unrelated -> 0
  prs = .comb2(ids, vec = TRUE)
  if(any(!mem$REL[prs]))
    return(debugReturn(0, debug = debug, indent = indent, comment = " (B3)"))

  # Sort: a > b > ..
  gp[[1]] = g = .mysort(g, decreasing = TRUE)

  # Lookup -> early return if previously computed
  s = paste(g, collapse = "-")
  if(!is.null(res <- mem$MEM[[s]])) {
    mem$ilook[L] = mem$ilook[L] + 1L
    return(debugReturn(res, debug = debug, indent = indent, comment = " (lookup)"))
  }

  ### Recurse
  mem$irec[L] = mem$irec[L] + 1L

  a = g[1]
  apar = a %% 10L # 1 (father) or 2 (mother)
  id1 = a %/% 10

  if(X && apar == 1L) {
    anc = mem$FIDX[id1]*10L + 2L
    kpNew = gipRepl1(gp, anc)
    res = psi1G(kpNew, X = X, mem = mem, debug = debug, indent = indent + 2)
  }
  else {
    anc = switch(apar, mem$FIDX[id1], mem$MIDX[id1]) * 10L + 1:2 # ancestral alleles
    kpNew1 = gipRepl1(gp, anc[1])
    kpNew2 = gipRepl1(gp, anc[2])
    res = (psi1G(kpNew1, X = X, mem = mem, debug = debug, indent = indent + 2) +
             psi1G(kpNew2, X = X, mem = mem, debug = debug, indent = indent + 2))/2
  }

  mem$MEM[[s]] = res
  debugReturn(res, debug = debug, indent = indent, comment = " (Computed)")
}


# Recursion function for the case of two groups with two gametes each.
# This is the only one needed for jacquard coeffs.
# Further generalised gamete probabilities are not currently implemented.
psi22 = function(gp, X = FALSE, mem = NULL, debug = FALSE, indent = 0) {
  if(debug)
    cat(strrep(" ", indent), gip2string(gp), "\n", sep = "")

  mem$i22 = mem$i22 + 1L
  g1 = gp[[1]]
  g2 = gp[[2]]

  ids1 = g1 %/% 10
  ids2 = g2 %/% 10

  # X males: pat -> mat
  if(X) {
    malePat1 = mem$SEX[ids1] == 1 & g1 %% 10 == 1
    if(any(malePat1)) g1[malePat1] = g1[malePat1] + 1L
    malePat2 = mem$SEX[ids2] == 1 & g2 %% 10 == 1
    if(any(malePat2)) g2[malePat2] = g2[malePat2] + 1L
  }

  if(g1[1] == g1[2] || g2[1] == g2[2]) {
    kpNew = if(g1[1] == g1[2]) gp[2] else gp[1]
    return(psi1G(kpNew, X = X, mem = mem, debug = debug, indent = indent + 2))
  }

  # Any pair unrelated -> 0
  if(!mem$REL[ids1[1], ids1[2]] || !mem$REL[ids2[1], ids2[2]])
    return(debugReturn(0, debug = debug, indent = indent, comment = " (B4)"))

  # (A:p, A:m) where A is not inbred -> 0
  if(ids1[1] == ids1[2] && mem$INB[ids1[1]] == 0 || ids2[1] == ids2[2] && mem$INB[ids2[1]] == 0)
    return(debugReturn(0, debug = debug, indent = indent, comment = " (B5)"))

  # All founders -> f1 * f2 (allows inbred founders)
  # NB: By now, the only case is (A:p, A:m), (B:p, B:m)
  if(all(mem$isFounder[c(ids1, ids2)])) {
    res = mem$INB[[ids1[1]]] * mem$INB[[ids2[1]]]
    return(debugReturn(res, debug = debug, indent = indent, comment = " (B6)"))
  }

  # Overlapping groups -> merge and recurse
  if(dup <- anyDuplicated.default(g <- c(g1, g2))) {
    kpNew = gp[1]
    kpNew[[1]] = g[-dup]
    return(psi1G(kpNew, X = X, mem = mem, debug = debug, indent = indent + 2))
  }

  # sort
  if(g1[1] < g1[2])
    gp[[1]] = g1 = g1[2:1]
  if(g2[1] < g2[2])
    gp[[2]] = g2 = g2[2:1]
  if(g1[1] < g2[1]) {
    tmp = g1; gp[[1]] = g1 = g2; gp[[2]] = g2 = tmp
  }

  # Lookup -> early return if previously computed
  s = paste(g1, g2, sep="-", collapse = "; ")
  if(!is.null(res <- mem$MEM[[s]])) {
    mem$ilook22 = mem$ilook22 + 1L
    return(debugReturn(res, debug = debug, indent = indent, comment = " (lookup)"))
  }

  # Recurse
  mem$irec22 = mem$irec22 + 1L

  a = g1[1]
  apar = a %% 10L # 1 (father) or 2 (mother)
  id1 = a %/% 10

  if(X && apar == 1L) {
    anc = mem$FIDX[id1]*10L + 2L
    kpNew = gipRepl1(gp, anc)
    res = psi22(kpNew, X = X, mem = mem, debug = debug, indent = indent + 2)
  }
  else {
    anc = switch(apar, mem$FIDX[id1], mem$MIDX[id1])*10L + 1:2
    kpNew1 = gipRepl1(gp, anc[1])
    kpNew2 = gipRepl1(gp, anc[2])
    res = (psi22(kpNew1, X = X, mem = mem, debug = debug, indent = indent + 2) +
             psi22(kpNew2, X = X, mem = mem, debug = debug, indent = indent + 2))/2
  }

  mem$MEM[[s]] = res
  debugReturn(res, debug = debug, indent = indent, comment = " (Computed)")
}
