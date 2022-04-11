twoLocusGeneralisedKinship = function(x, locus1, locus2, rho, verbose = FALSE, debug = FALSE) {

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  x = foundersFirst(x)

  # Setup memoisation
  mem = initialiseTwoLocusMemo(x, rho, counters = c("i", "itriv", "iimp", "ifound", "ilook", "irec"))

  kc = kin2L(x, locus1, locus2)

  res = genKin2L(kc, mem, indent = ifelse(debug, 0, NA))

  if(verbose)
    printCounts2(mem)

  res
}


genKin2L = function(kin, mem, indent = 0) {
  mem$i = mem$i + 1

  print.kin2L(kin, indent = indent)
  kin = kinReduce(kin)

  ll = lengths(kin)
  if(sum(ll) == 0) {
    mem$itriv = mem$itriv + 1
    return(printAndReturn(1, indent, comment = " (empty)"))
  }

  if(all(ll <= 1)) {
    triv1 = ll[1] == 0 || length(kin$locus1[[1]]$from) == 1
    triv2 = ll[2] == 0 || length(kin$locus2[[1]]$from) == 1
    if(triv1 && triv2) {
      mem$itriv = mem$itriv + 1
      return(printAndReturn(1, indent, comment = " (trivial)"))
    }
  }

  if(kinImpossible(kin, mem))
    return(printAndReturn(0, indent, comment = " (impossible)"))

  if(boundary_test(kin, mem)) {
    val = boundary_value(kin, mem)
    return(printAndReturn(val, indent, comment = " (boundary 3)"))
  }

  kin = sort_kin2L(kin)

  # Lookup in array; compute if necessary.
  kinStr = toString.default(kin)
  val = mem$k2[[kinStr]]

  if(!is.null(val)) {
    mem$ilook = mem$ilook + 1
    mem$k2[[kinStr]] = val
    return(printAndReturn(val, indent, comment = " (lookup)"))
  }

  mem$irec = mem$irec + 1
  rho = mem$rho

  # Pivot individual and parents
  a = kin$locus1[[1]]$from[1]
  f = mem$FIDX[a]
  m = mem$MIDX[a]

  # Targets of pivot: list of 4 vectors (loc1.g1, loc1.g2, loc2.g1, loc2.g2)
  pivTargets = lapply(c(kin$locus1[1:2], kin$locus2[1:2]),
                      function(g) g$to[g$from == a])

  r = length(pivTargets[[1]])
  s = length(pivTargets[[2]])
  t = length(pivTargets[[3]])
  u = length(pivTargets[[4]])


  un = length(unique.default(unlist(pivTargets, use.names = FALSE)))
  ev = length(c(intersect(pivTargets[[1]], pivTargets[[3]]),
           intersect(pivTargets[[2]], pivTargets[[4]])))
  od = length(c(intersect(pivTargets[[1]], pivTargets[[4]]),
         intersect(pivTargets[[2]], pivTargets[[3]])))

  if(!is.na(indent)) {
    message(sprintf("%sRecurse: a = %d; father = %d; mother = %d",
                    strrep(" ", indent), a,f,m))
  }
  # Wrapper of genKin2L to save typing in the recursions
  recu = function(k) genKin2L(k, mem, indent = indent +2)

  # Meiosis indicators used in recursions (separates f and m in case of selfing)
  a.f = 100*a + 1
  a.m = 100*a + 2

  # Recurse!
  if(s + t + u == 0) {
    A1 = .5^r
    res =
      A1 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f))) +
      A1 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m)))

    if(r > 1) {
      B1 = (1 - 2 * (.5)^r)
      res = res +
        B1 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f, a.m))))
    }
  }
  else if(t + u == 0) {
    A2 = .5^(r + s)
    res =
      A2 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m))) +
      A2 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f)))
  }
  else if(s + u == 0) {
    A3 = .5^un * (1-rho)^ev
    B3 = .5^un * rho^ev
    res =
      A3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f), loc2Rep = list(from1 = f, to1 = a.f))) +
      B3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f), loc2Rep = list(from1 = m, to1 = a.m))) +
      B3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m), loc2Rep = list(from1 = f, to1 = a.f))) +
      A3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m), loc2Rep = list(from1 = m, to1 = a.m)))

    C3 = D3 = E3 = 0
    R.ev = (1-rho)^ev + rho^ev

    if(r > 1) {
      C3 = .5^t - .5^un * R.ev
      res = res +
        C3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)), loc2Rep = list(from1 = f, to1 = a.f))) +
        C3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)), loc2Rep = list(from1 = m, to1 = a.m)))
    }
    if(t > 1) {
      D3 = .5^r - .5^un * R.ev
      res = res +
        D3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)))) +
        D3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m))))
    }
    if(r > 1 && t > 1) {
      E3 = 1 - 2*.5^t - 2 * .5^r + 2*.5^un*R.ev
      res = res +
        E3 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m))))
    }

    # Check that coeffs sum to 1
    if((SM <- 2*A3 + 2*B3 + 2*C3 + 2*D3 + E3) != 1)
      stop2("Something wrong in case 's = u = 0; coefs sum to: ", SM)
  }
  else if(u == 0) {
    A4 = .5^un * (1-rho)^ev * rho^od
    B4 = .5^un * (1-rho)^od * rho^ev
    res =
      A4 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = f, to1 = a.f))) +
      B4 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = f, to1 = a.f))) +
      B4 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = m, to1 = a.m))) +
      A4 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = m, to1 = a.m)))
    if(t > 1) {
      C4 = .5^(r+s) - .5^un * ((1-rho)^ev * rho^od + (1-rho)^od * rho^ev)
      res = res +
        C4 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)))) +
        C4 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m))))
    }
  }
  else {
    A5 = .5^un * (1-rho)^ev * rho^od
    B5 = .5^un * (1-rho)^od * rho^ev
    res =
      A5 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m))) +
      B5 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f))) +
      B5 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m))) +
      A5 * recu(kinReplace(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f)))
  }

  mem$k2[[kinStr]] = res
  return(printAndReturn(res, indent))
}


printCounts2 = function(mem) {
  tot_time = format(Sys.time()-mem$st, digits = 3)

  msg = glue::glue("
                   Function calls: {mem$i}
                     Trivial:      {mem$itriv}
                     Impossible:   {mem$iimp}
                     Boundary:     {mem$ifound}
                     Lookups:      {mem$ilook}
                     Recursions:   {mem$irec}
                   Total time: {tot_time}")
  message(msg)
}


