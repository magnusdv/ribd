# Boundary condition 0: Any identical from>to in different groups?
boundary0_test = function(kin, mem) {
  for(i in 1:2) {
    meioses = lapply(kin[[i]], function(g) paste(g$from, g$to, sep=">"))
    if(anyDuplicated(unlist(meioses))) {
      mem$b0 = mem$b0 + 1
      return(T)
    }
  }
  return(F)
}

# Boundary condition 1: Anyone in >2 groups?
boundary1_test = function(kin, mem) {
  for(i in 1:2) {
    uniq.ids = lapply(kin[[i]], function(g) unique.default(g$from))
    if(length(uniq.ids) == 0)
      next
    tab = tabulate(unlist(uniq.ids))
    if(any(tab > 2)) {
      mem$b1 = mem$b1 + 1
      return(T)
    }
  }
  return(F)
}


# Boundary condition 2: Any group with 2 unrelated indivs?
boundary2_test = function(kin, mem) {
  k1 = mem$k1
  for(g in c(kin[[1]], kin[[2]])) {
    ids = unique.default(g$from)
    if(length(ids) < 2)
      next
    pairs_mat = comb2(ids, vec = T)
    if(any(k1[pairs_mat] == 0)) {
      mem$b2 = mem$b2 + 1
      return(T)
    }
  }
  return(F)
}

boundary3_test = function(kin, mem) {
  loc1 = unlist(lapply(kin$locus1, function(g) g$from))
  loc2 = unlist(lapply(kin$locus2, function(g) g$from))
  if(all(mem$isFounder[c(loc1, loc2)])) {
    mem$b3 = mem$b3 + 1
    return(T)
  }
  return(F)
}

# If pivot and everybody else is founder,
# and neither boundary 1 or 2
boundary3_value = function(kin, mem) {

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

    # Only present in locus1
    if(length(L2a) == 0) {
      res = res * 0.5 ^ (length(unlist(L1a))/2 - 1) # dividing by 2: counting both from/to
      next
    }

    # Only present in locus2
    if(length(L1a) == 0) {
      res = res * 0.5 ^ (length(unlist(L2a))/2 - 1)
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
  nAll = length(unique.default(unlist(x, use.names = F)))

  loc1.g1 = x[[1]]
  loc1.g2 = x[[2]]
  loc2.g1 = x[[3]]
  loc2.g2 = x[[4]]
  even = c(intersect(loc1.g1, loc2.g1), intersect(loc1.g2, loc2.g2))
  odd  = c(intersect(loc1.g1, loc2.g2), intersect(loc1.g2, loc2.g1))

  list(n = nAll, even = length(even), odd = length(odd))
}
