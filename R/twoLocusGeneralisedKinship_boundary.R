
kinImpossible = function(kin, mem) {

  # Loop over the two loci
  for(i in 1:2) {

    locusGroups = kin[[i]]
    if(length(locusGroups) == 0)
      next

    # List unique source indivs in each group
    uniq_sources = lapply(locusGroups, function(g) unique.default(g$from))

    # Boundary condition: Any group with 2 unrelated indivs?
    k1 = mem$k1
    for(s in uniq_sources[lengths(uniq_sources) > 1]) {
      pairs_mat = .comb2(s, vec = TRUE)
      if(any(k1[pairs_mat] == 0)) {
        mem$iimp = mem$iimp + 1
        return(TRUE)
      }
    }

    # Boundary condition: Anyone in >2 groups?
    tab = tabulate(unlist(uniq_sources, use.names = FALSE))
    if(any(tab > 2)) {
      mem$iimp = mem$iimp + 1
      return(TRUE)
    }

    # Boundary condition 0: Any identical from>to in different groups?
    meioses = lapply(locusGroups, function(g) 1000*g$to + g$from)
    if(anyDuplicated.default(unlist(meioses, use.names = FALSE))) {
      mem$iimp = mem$iimp + 1
      return(TRUE)
    }
  }

  # If none of the boundary conditions are met, return FALSE
  return(FALSE)
}


# Check: Are all sources founders?
boundary_test = function(kin, mem) {
  isFou = mem$isFounder

  loc1 = unlist(lapply(kin$locus1, function(g) g$from), use.names = FALSE)
  if(!all(isFou[loc1]))
    return(FALSE)

  loc2 = unlist(lapply(kin$locus2, function(g) g$from), use.names = FALSE)
  if(!all(isFou[loc2]))
    return(FALSE)

  mem$ifound = mem$ifound + 1
  return(TRUE)
}

# Return value in boundary case
# Assuming `kinImpossible(kin)` is FALSE!
boundary_value = function(kin, mem) {

  L1 = kin$locus1
  L2 = kin$locus2

  # Each group consists of a single founder
  # (possibly repeated, but then with different targets)
  ids1 = vapply(L1, function(g) g$from[1], 1)
  ids2 = vapply(L2, function(g) g$from[1], 1)

    # Initialise result
  res = 1

  # Loop over all involved founders
  for(a in unique.default(c(ids1, ids2))) {
    L1a = L1[ids1 == a]
    L2a = L2[ids2 == a]

    # Completely inbred?
    if(mem$isCompletelyInbred[a]) {
      if(length(L1a) > 1 || length(L2a) > 1) { # Impossible!
        res = 0
        break
      }
      res = res * 1  # Otherwise trivial
      next
    }

    ### By now a can be assumed to be outbred ###

    # Only present in locus1
    if(length(L2a) == 0) {
      res = res * 0.5 ^ (length(unlist(L1a, use.names = FALSE))/2 - 1) # dividing by 2: counting both from/to
      next
    }

    # Only present in locus2
    if(length(L1a) == 0) {
      res = res * 0.5 ^ (length(unlist(L2a, use.names = FALSE))/2 - 1)
      next
    }

    ### In both loci
    r = mem$rho

    # Extract targets from groups 1 and 2 (poss empty) at each locus
    to1 = lapply(L1a, function(g) g$to)
    if(length(to1) == 1) to1[[2]] = integer(0)

    to2 = lapply(L2a, function(g) g$to)
    if(length(to2) == 1) to2[[2]] = integer(0)

    # All unique targets except 1 contribute with 1/2.
    # In addition, those in both loci get a factor [r^n * (1-r)^m + r^m * (1-r)^n]

    # Exponents m and n are determined by parity considerations
    parit = parityCounts(c(to1, to2))
    even = parit$even
    odd = parit$odd
    n = parit$n
    #cat("even: ", even, "odd: ", odd, "onelocus: ", oneLocus, "\n")
    #cat(sprintf("0.5 ^(%d) * (r^%d * (1-r)^%d + (1-r)^%d * r^%d)\n", n - 1, even, odd, even, odd))
    res = res * 0.5 ^(n - 1) * (r^even * (1-r)^odd + (1-r)^even * r^odd)
  }

  res
}

# Input: list of 4 vectors
parityCounts = function(x) {
  nAll = length(unique.default(unlist(x, use.names = FALSE)))

  loc1.g1 = x[[1]]
  loc1.g2 = x[[2]]
  loc2.g1 = x[[3]]
  loc2.g2 = x[[4]]

  # TODO: Presumably this never happens after kinReduce.
  if(any(sapply(x, anyDuplicated.default) > 0))
    stop2("Duplicated targets! Please contact maintainer")

  even = c(.myintersect(loc1.g1, loc2.g1), .myintersect(loc1.g2, loc2.g2))
  odd  = c(.myintersect(loc1.g1, loc2.g2), .myintersect(loc1.g2, loc2.g1))

  list(n = nAll, even = length(even), odd = length(odd))
}
