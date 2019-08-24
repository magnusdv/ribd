twoLocusGeneralisedKinship = function(x, locus1, locus2, rho, verbose = F, debug = F) {

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  x = foundersFirst(x)

  # Setup memoisation
  mem = initialiseTwoLocusMemo(x, rho, counters = c("i", "ilook", "ir", "b0", "b1", "b2", "b3"))

  kc = kin2L(x, locus1, locus2)

  res = genKin2L(kc, mem, indent = ifelse(debug, 0, NA))

  if(verbose)
    printCounts2(mem)

  res
}


genKin2L = function(kin, mem, indent = 0) {
  mem$i = mem$i + 1

  print(kin, indent = indent)
  kin = kinReduce(kin)

  ll = lengths(kin)
  if(sum(ll) == 0)
    return(printAndReturn(1, indent, comment = " (empty)"))

  if(all(ll <= 1)) {
    triv1 = ll[1] == 0 || length(kin$locus1[[1]]$from) == 1
    triv2 = ll[2] == 0 || length(kin$locus2[[1]]$from) == 1
    if(triv1 && triv2)
      return(printAndReturn(1, indent, comment = " (trivial)"))
  }

  if(boundary0_test(kin, mem))
    return(printAndReturn(0, indent, comment = " (boundary 0)"))

  if(boundary1_test(kin, mem))
    return(printAndReturn(0, indent, comment = " (boundary 1)"))

  if(boundary2_test(kin, mem))
    return(printAndReturn(0, indent, comment = " (boundary 2)"))

  if(boundary3_test(kin, mem)) {
    val = boundary3_value(kin, mem)
    return(printAndReturn(val, indent, comment = " (boundary 3)"))
  }

  kin = sort_kin2L(kin)

  # Lookup in array; compute if necessary.
  kinStr = toString(kin)
  val = mem$k2[[kinStr]]

  if(!is.null(val)) {
    mem$ilook = mem$ilook + 1
    mem$k2[[kinStr]] = val
    return(printAndReturn(val, indent, comment = " (lookup)"))
  }

  mem$ir = mem$ir + 1
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


  un = length(unique.default(unlist(pivTargets, use.names = F)))
  ev = length(c(intersect(pivTargets[[1]], pivTargets[[3]]),
           intersect(pivTargets[[2]], pivTargets[[4]])))
  od = length(c(intersect(pivTargets[[1]], pivTargets[[4]]),
         intersect(pivTargets[[2]], pivTargets[[3]])))

  if(!is.na(indent)) {
    cat(sprintf("%sr,s,t,u = %d,%d,%d,%d; uniq = %d; even = %d; odd = %d\n", strrep(" ", indent), r,s,t,u,un,ev,od))
    message(sprintf("%sRecurse: a = %d; father = %d; mother = %d", strrep(" ", indent), a,f,m))
  }
  # Wrapper of genKin2L to save typing in the recursions
  recu = function(k) genKin2L(k, mem, indent = indent +2)

  # Meiosis indicators used in recursions (separates f and m in case of selfing)
  a.f = 100*a + 1
  a.m = 100*a + 2

  # Recurse!
  if(s + t + u == 0) {
    res =
      .5^r * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f))) +
      .5^r * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m)))

    if(r > 1)
      res = res +
        (1 - 2 * (.5)^r) * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f, a.m))))

    return(printAndReturn(res, indent))
  }

  if(t + u == 0) {
    A.2 = .5^(r + s)
    res =
      A.2 * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m))) +
      A.2 * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f)))

    return(printAndReturn(res, indent))

  }
  if(s + u == 0) {
    A1 = .5^un * (1-rho)^ev
    A2 = .5^un * rho^ev
    res =
      A1 * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f), loc2Rep = list(from1 = f, to1 = a.f))) +
      A2 * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f), loc2Rep = list(from1 = m, to1 = a.m))) +
      A2 * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m), loc2Rep = list(from1 = f, to1 = a.f))) +
      A1 * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m), loc2Rep = list(from1 = m, to1 = a.m)))

    B = C = D = 0
    R.ev = (1-rho)^ev + rho^ev

    if(r > 1) {
      B = .5^t - .5^un * R.ev
      res = res +
        B * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)), loc2Rep = list(from1 = f, to1 = a.f))) +
        B * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)), loc2Rep = list(from1 = m, to1 = a.m)))
    }
    if(t > 1) {
      C = .5^r - .5^un * R.ev
      res = res +
        C * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)))) +
        C * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m))))
    }
    if(r > 1 && t > 1) {
      D = 1 - 2*.5^t - 2 * .5^r + 2*.5^un*R.ev
      res = res +
        D * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m))))
    }

    # Check that coeffs sum to 1
    if((SM <- 2*A1 + 2*A2 + 2*B + 2*C + D) != 1)
      stop2("Something wrong in case 's=u=0; coefs sum to: ", SM)
    return(printAndReturn(res, indent))
  }

  if(u == 0) {
    A = .5^un * (1-rho)^ev * rho^od
    B = .5^un * (1-rho)^od * rho^ev
    res =
      A * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = f, to1 = a.f))) +
      B * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = f, to1 = a.f))) +
      B * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = m, to1 = a.m))) +
      A * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = m, to1 = a.m)))
    if(t > 1) {
      C = .5^(r+s) - .5^un * ((1-rho)^ev * rho^od + (1-rho)^od * rho^ev)
      res = res +
        C * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m)))) +
        C * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = c(f,m), to1 = c(a.f,a.m))))
    }

    return(printAndReturn(res, indent))
  }

  A = .5^un * (1-rho)^ev * rho^od
  B = .5^un * (1-rho)^od * rho^ev
  res =
    A * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m))) +
    B * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m), loc2Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f))) +
    B * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = f, to1 = a.f, from2 = m, to2 = a.m))) +
    A * recu(kinRepl(kin, id = a, loc1Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f), loc2Rep = list(from1 = m, to1 = a.m, from2 = f, to2 = a.f)))

  return(printAndReturn(res, indent))
}


printCounts2 = function(mem) {
  initsecs = sprintf("%.2f", mem$initTime)
  totsecs = sprintf("%.1f", Sys.time()-mem$st)

  msg = glue::glue("
                   Function calls: {mem$i}
                   Lookups:        {mem$ilook}
                   Recursions:     {mem$ir}
                   Boundary cases
                     B0: {mem$b0}
                     B1: {mem$b1}
                     B2: {mem$b2}
                     B3: {mem$b3}
                   Initialization: {initsecs} seconds
                   Total time:     {totsecs} seconds")
  print(msg)
}


