#####################################################################
# Recursive algorithm for identity coefficients by Lange & Sinsheimer
#
# Note: Different algorithms for condensed and detailed coeffs
#####################################################################


identity_LS = function(x, ids, Xchrom = FALSE, detailed = FALSE, self = FALSE, mem = NULL, verbose = FALSE) {

  x = prepPed(x, addpar = ids, Xchrom = Xchrom)
  pairs = prepIds(x, ids, self = self)

  if(is.null(mem))
    mem = memoIdentity(x, method = "LS")

  # Recursion wrapper to simplify typing
  recu = function(..., debug = FALSE) {
    recurse_LS(kinpat(x, list(...)), X = Xchrom, mem = mem, debug = debug)
  }

  # Detailed identity states (following Figure 6.2 (page 105), Jacquard 1974)
  PhiLS = function(id1, id2) {
    c(s1 = recu(c(p=id1, m=id1, p=id2, m=id2)),
      s2 = recu(c(p=id1, m=id1, p=id2), c(m=id2)),
      s3 = recu(c(p=id1, m=id1, m=id2), c(p=id2)),
      s4 = recu(c(p=id1, p=id2, m=id2), c(m=id1)),
      s5 = recu(c(m=id1, p=id2, m=id2), c(p=id1)),
      s6 = recu(c(p=id1, m=id1), c(p=id2, m=id2)),
      s7 = recu(c(p=id1, m=id1), c(p=id2), c(m=id2)),
      s8 = recu(c(p=id1), c(m=id1), c(p=id2, m=id2)),
      s9 = recu(c(p=id1, p=id2), c(m=id1, m=id2)),
      s10 = recu(c(p=id1, p=id2), c(m=id1), c(m=id2)),
      s11 = recu(c(p=id1), c(p=id2), c(m=id1, m=id2)),
      s12 = recu(c(p=id1, m=id2), c(m=id1, p=id2)),
      s13 = recu(c(p=id1, m=id2), c(m=id1), c(p=id2)),
      s14 = recu(c(p=id1), c(m=id2), c(m=id1, p=id2)),
      s15 = recu(c(p=id1), c(m=id1), c(m=id2), c(p=id2)))
  }

  # Compute Phi for each pair
  j = vapply(pairs, function(p) PhiLS(p[1], p[2]), FUN.VALUE = numeric(15))

  if(verbose)
    printMemInfo(mem)

  # Data frame
  res = jmat2df(t.default(j), pairs)

  if(!detailed)
    res = detailed2condensed(res)

  res
}

gKinship_LS = function(x, kp, Xchrom = FALSE, mem, debug = FALSE) {
  recurse_LS(kp, X = Xchrom, mem = mem, debug = debug)
}

# Recursion method of Lange & Sinsheimer
recurse_LS = function(kp, X = FALSE, mem = NULL, debug = FALSE, indent = 0) {
  mem$i = mem$i + 1
  kp = kinpatReduce(kp, deterministic = TRUE)

  # If no longer deterministic, switch to random algorithm (WL)
  if(!isDeterministic(kp))
    return(recurse_WL(kp, X = X, mem, debug = debug, indent = indent))

  if(debug)
    cat(strrep(" ", indent), kinpat2string(kp), "\n", sep = "")

  # B0: Trivial pattern
  if(sum(lengths(kp)) <= 1)  {
    mem$B0 = mem$B0 + 1
    return(debugReturn(1, debug = debug, indent = indent, comment = " (B0)"))
  }

  # B2: Any group with 2 unrelated indivs?
  uniqList = lapply(kp, function(g) unique.default(g %/% 10))
  for(s in uniqList) {
    if(length(s) < 2)
      next
    intrapairs = comb2(s, vec = TRUE)
    if(any(!mem$REL[intrapairs])) { # any pair not related?
      mem$B2 = mem$B2 + 1
      return(debugReturn(0, debug = debug, indent = indent, comment = " (B2)"))
    }
  }

  # B1: Anyone in >2 groups
  uniqVec = unlist(uniqList, use.names = FALSE)
  tab = tabulate(uniqVec)
  if(any(tab > 2)) {
    mem$B1 = mem$B1 + 1
    return(debugReturn(0, debug = debug, indent = indent, comment = " (B1)"))
  }

  # B3: Only founders left: This should never happen, since original founders have added parents

  # Wrapper to save typing in the recursions
  recu = function(k) recurse_LS(k, X = X, mem = mem, debug = debug, indent = indent +2)

  # If 2 (for simplicity) unrelated groups: Factorise
  if(length(kp) == 2) {
    u1 = uniqList[[1]]; u2 = uniqList[[2]]
    interpairs = cbind(rep(u1, each = length(u2)), rep(u2, length(u1)))
    if(!any(mem$REL[interpairs])) { # any pair not related?
      res = recu(`[[<-`(kp, 2, NULL)) * recu(`[[<-`(kp, 1, NULL))
      return(debugReturn(res, debug = debug, indent = indent, comment = " Factorised"))
    }
  }

  ### Additional restrictions for deterministic patterns

  # Restriction D2: Anyone with m&p in same group AND (in multiple groups OR not inbred)
  if(restrictionD2(kp, tab = tab, inb = mem$INB)) {
    mem$D2 = mem$D2 + 1
    return(debugReturn(0, debug = debug, indent = indent, comment = " (D2)"))
  }

  # Restriction D1: Anyone with paternal (resp. maternal) allele in multiple blocks?
  if(any(tab > 1) && restrictionD1(kp)) {
    mem$D1 = mem$D1 + 1
    return(debugReturn(0, debug = debug, indent = indent, comment = " (D1)"))
  }

  kp = kinpatSort(kp)

  # Lookup in array; compute if necessary.
  kinStr = paste(kp, collapse = ", ")
  val = mem$MEM2[[kinStr]]

  if(!is.null(val)) {
    mem$ilook = mem$ilook + 1
    return(debugReturn(val, debug = debug, indent = indent, comment = " (lookup)"))
  }

  mem$irec = mem$irec + 1

  g1 = kp[[1]]
  g2 = if(length(kp) > 1) kp[[2]] else NULL

  ids1 = g1 %/% 10
  ids2 = g2 %/% 10

  # Pivot individual (after sorting: highest)
  pivot = ids1[1]
  fa = mem$FIDX[pivot]
  mo = mem$MIDX[pivot]
  famo = c(fa, mo)

  # Number of times the pivot occurs in first block (s) and second block (t)
  s = sum(ids1 == pivot)
  t = sum(ids2 == pivot)

  # Deterministic sampling of pivot
  par1 = g1[ids1 == pivot] %% 10
  det1 = par1[par1 > 0]  # vector of length 0, 1 or 2
  ndet1 = length(det1)
  s = s - ndet1

  # Recurrence rule 1a,b,c
  if(t == 0) {
    A1 = .5^s
    res = switch(ndet1 + 1,
                 # Rule 1a
                 A1 * recu(kinpatReplaceDet(kp, id = pivot, rep1 = fa)) +
                   A1 * recu(kinpatReplaceDet(kp, id = pivot, rep1 = mo)) +
                   if(s > 1) (1 - 2*A1) * recu(kinpatReplaceDet(kp, id = pivot, rep1 = famo)) else 0,

                 # Rule 1b
                 A1 * recu(kinpatReplaceDet(kp, id = pivot, rep1 = famo[det1])) +
                   if(s > 0) (1 - A1) * recu(kinpatReplaceDet(kp, id = pivot, rep1 = famo)) else 0,

                 # Rule 1c
                 res = recu(kinpatReplaceDet(kp, id = pivot, rep1 = famo))
    )
  }
  else {  # t > 0

    # Deterministic sampling in second group
    par2 = g2[ids2 == pivot] %% 10
    det2 = par2[par2 > 0]  # vector of length 0, 1 or 2
    ndet2 = length(det2)
    t = t - ndet2

    A2 = .5^(s + t)

    if(ndet1 + ndet2 == 0) { # Rule 2a
      res =
        A2 * recu(kinpatReplaceDet(kp, id = pivot, rep1 = fa, rep2 = mo)) +
        A2 * recu(kinpatReplaceDet(kp, id = pivot, rep1 = mo, rep2 = fa))
    }
    else {  # Rules 2b and 2c
      # NB: This must also work when fa = mo (selfing!)
      reps = if(ndet1 == 1 && det1 == 1 || ndet2 == 1 && det2 == 2) famo else rev.default(famo)
      res = A2 * recu(kinpatReplaceDet(kp, id = pivot, rep1 = reps[1], rep2 = reps[2]))
    }
  }

  mem$MEM2[[kinStr]] = res
  return(debugReturn(res, debug = debug, indent = indent, comment = " (Computed)"))
}


# Restriction D1: Anyone with paternal (resp. maternal) allele in multiple blocks?
restrictionD1 = function(kp) {
  # Assumes that kp is reduced, i.e., nobody has repeated det-sampling *within* a block
  vec = unlist(kp, use.names = FALSE)
  anyDuplicated.default(vec[vec %% 10 > 0])
}

# Restriction D2: Anyone with m&p in same group AND (in multiple groups OR not inbred)
restrictionD2 = function(kp, tab, inb) {
  for(g in kp) {
    gdet = g[g %% 10 > 0]
    if(length(gdet) < 2) next
    idsdet = gdet %/% 10
    mp = idsdet[duplicated.default(idsdet)]
    if(any(tab[mp] > 1) || any(inb[mp] == 0))
      return(TRUE)
  }
  FALSE
}

