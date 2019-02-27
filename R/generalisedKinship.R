gKinship3 = function(x, ids, sparse = NA, verbose = TRUE) {

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, verbose = verbose)

  # Compute
  res = phi3(ids_int[1], ids_int[2], ids_int[3], mem)

  if(verbose)
    printCounts(mem)

  res
}


gKinship4 = function(x, ids, sparse = NA, verbose = TRUE) {

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, verbose = verbose)

  res = phi4(ids_int[1], ids_int[2], ids_int[3], ids_int[4], mem)
  if(verbose)
    printCounts(mem)

  res
}

gKinship22 = function(ped, ids, sparse = NA, verbose = TRUE) {
  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, verbose = verbose)

  res = phi22(ids_int[1], ids_int[2], ids_int[3], ids_int[4], mem)
  if(verbose)
    printCounts(mem)

  res
}

phi2 = function(a, b, mem) {
  mem$i2 = mem$i2 + 1
  if(a*b == 0) return(0)

  mem$KIN2[[a, b]]
}

phi3 = function(a, b, c, mem) {

  mem$i3 = mem$i3 + 1
  if(a*b*c == 0) return(0)

  ANC = mem$ANC
  FIDX = mem$FIDX
  MIDX = mem$MIDX

  if(!all(ANC[a,b],ANC[b,c],ANC[a,c])) return(0)

  # Sort: a>=b>=c
  if(a < b) {tmp=a; a=b; b=tmp}
  if(b < c) {tmp=b; b=c; c=tmp}
  if(a < b) {tmp=a; a=b; b=tmp}

  # Recursion function (called only if necessary)
  # Assumes a,b,c sorted
  phi3_recurse = function(a, b, c) {
    mem$i3r = mem$i3r + 1

    if(a == b && a == c)
      return((1 + 3*phi2(FIDX[a], MIDX[a], mem = mem))/4)

    if(a == b)
      return((phi2(a, c, mem = mem) + phi3(FIDX[a], MIDX[a], c, mem = mem))/2)

    return((phi3(FIDX[a], b, c, mem = mem) + phi3(MIDX[a], b, c, mem = mem))/2)
  }

  # Lookup in array; compute if necessary.
  res = mem$KIN3[[a,b,c]]
  if(is.na(res))
    res = mem$KIN3[[a,b,c]] = phi3_recurse(a, b, c)

  res
}

phi4 = function(a, b, c, d, mem) {

  mem$i4 = mem$i4 + 1
  if(a*b*c*d == 0) return(0)

  ANC = mem$ANC
  FIDX = mem$FIDX
  MIDX = mem$MIDX
  if(!all(ANC[a,b],ANC[b,c],ANC[c,d],ANC[b,d],ANC[a,c],ANC[a,d])) return(0)

  # Sort: a >= b >= c >= d
  if(a < b) {tmp=a; a=b; b=tmp}
  if(a < c) {tmp=a; a=c; c=tmp}
  if(a < d) {tmp=a; a=d; d=tmp}
  if(b < c) {tmp=b; b=c; c=tmp}
  if(b < d) {tmp=b; b=d; d=tmp}
  if(c < d) {tmp=c; c=d; d=tmp}

  # Recursion function (called only if necessary)
  # Assumes a,b,c,d sorted
  phi4_recurse = function(a, b, c, d) {
    mem$i4r = mem$i4r + 1

    if(a == b && a == c && a == d)
      return((1 + 7*phi2(FIDX[a], MIDX[a], mem = mem))/8)

    if(a == b && a == c)
      return((phi2(a, d, mem = mem) + 3*phi3(FIDX[a], MIDX[a], d, mem = mem))/4)

    if(a == b)
      return((phi3(a, c, d, mem = mem) + phi4(FIDX[a], MIDX[a], c, d, mem = mem))/2)

    return((phi4(FIDX[a], b, c, d, mem = mem) + phi4(MIDX[a], b, c, d, mem = mem))/2)
  }

  # Lookup in array; compute if necessary.
  res = mem$KIN4[[a,b,c,d]]
  if(is.na(res))
    res = mem$KIN4[[a,b,c,d]] = phi4_recurse(a, b, c, d)

  res
}

phi22 = function(a, b, c, d, mem = NULL) {
  mem$i22 = mem$i22 + 1
  if(a*b*c*d == 0) return(0)

  ANC = mem$ANC
  FIDX = mem$FIDX
  MIDX = mem$MIDX
  if(!(ANC[a,b] && ANC[c,d])) return(0)

  # Sort: a >= b,c,d; c >= d; if(a==c) then b>=d
  s = c(min(a,b), max(a,b), min(c,d), max(c,d)) # d,c,b,a
  if(s[4] < s[2] || (s[4] == s[2] && s[3] < s[1]))
    s[] = s[c(3,4,1,2)]
  a = s[4]; b = s[3]; c = s[2]; d = s[1]

  # Recursion function (called only if necessary)
  # Assumes a,b,c sorted
  phi22_recurse = function(a, b, c, d) {
    mem$i22r = mem$i22r + 1

    if(a == b && a == c && a == d)
      return((1 + 3*phi2(FIDX[a], MIDX[a], mem = mem))/4)

    if(a == b && a == c)
      return((phi2(a, d, mem = mem) + phi3(FIDX[a], MIDX[a], d, mem = mem))/2)

    if(a == b) { #NB modification to allow inbred founders!
      if(mem$isFounder[a])
        return(0.5*phi2(c, d, mem = mem)*(1 + mem$founderInb[a]))

      return((phi2(c, d, mem = mem) + phi22(FIDX[a], MIDX[a], c, d, mem = mem))/2)
    }

    if(a == c)
      return((2*phi3(a, b, d, mem = mem) +
                phi22(FIDX[a], b, MIDX[a], d, mem = mem) +
                phi22(MIDX[a], b, FIDX[a], d, mem = mem))/4)

    return((phi22(FIDX[a], b, c, d, mem = mem) + phi22(MIDX[a], b, c, d, mem = mem))/2)
  }

  # Lookup in array; compute if necessary.
  res = mem$KIN22[[a,b,c,d]]
  if(is.na(res))
    res = mem$KIN22[[a,b,c,d]] = phi22_recurse(a, b, c, d)

  res
}
