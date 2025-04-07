
identity_Karigl = function(x, ids, Xchrom = FALSE, sparse = 30, self = FALSE, verbose = FALSE) {

  x = prepPed(x)
  pairs = prepIds(x, ids, self = self, int = TRUE)

  # Setup memoisation
  mem = memoIdentity(x, Xchrom = Xchrom, method = "K", sparse = sparse, maxId = max(unlist(pairs)),
                     counters = c("i2", "i3", "i4", "i22", "i2r", "i3r", "i4r", "i22r"), verbose = verbose)

  # M9 = matrix(c(
  #   1,1,1,1,1,1,1,1,1,
  #   2,2,2,2,1,1,1,1,1,
  #   2,2,1,1,2,2,1,1,1,
  #   4,0,2,0,2,0,2,1,0,
  #   8,0,4,0,2,0,2,1,0,
  #   8,0,2,0,4,0,2,1,0,
  #   16,0,4,0,4,0,2,1,0,
  #   4,4,2,2,2,2,1,1,1,
  #   16,0,4,0,4,0,4,1,0), byrow = TRUE, ncol = 9)

  Minv = matrix(ncol = 9, byrow = TRUE, c(
    0, 0, 0, 0.25,-0.25,-0.25, 0.25, 0,   0,
    1,-1,-1,-0.25, 0.25, 0.25,-0.25, 1,   0,
    0, 0, 0,   -1,    1,  0.5, -0.5, 0,   0,
   -2, 2, 1,    1,   -1, -0.5,  0.5,-1,   0,
    0, 0, 0,   -1,  0.5,    1, -0.5, 0,   0,
   -2, 1, 2,    1, -0.5,   -1,  0.5,-1,   0,
    0, 0, 0,    0,    0,    0, -0.5, 0, 0.5,
    0, 0, 0,    4,   -2,   -2,    2, 0,-1.0,
    4,-2,-2,   -4,    2,    2, -1.5, 1, 0.5
  ))

  RHS = vapply(pairs, function(p) {
    id1 = p[1]; id2 = p[2]
    c(1,
      2 * phi2(id1, id1, X = Xchrom, mem = mem),
      2 * phi2(id2, id2, X = Xchrom, mem = mem),
      4 * phi2(id1, id2, X = Xchrom, mem = mem),
      8 * phi3(id1, id1, id2, X = Xchrom, mem = mem),
      8 * phi3(id1, id2, id2, X = Xchrom, mem = mem),
      16 * phi4(id1, id1, id2, id2, X = Xchrom, mem = mem),
      4 * phi22(id1, id1, id2, id2, X = Xchrom, mem = mem),
      16 * phi22(id1, id2, id1, id2, X = Xchrom, mem = mem))
  }, FUN.VALUE = numeric(9))

  # Compute identity coefficients: matrix with 9 rows
  j = Minv %*% RHS

  if(verbose)
    printMemInfo(mem)

  # Build result data frame
  res = jmat2df(t.default(j), pairs, labs = labels(x))

  if(Xchrom)
    res = Xmask(x, res)

  res
}


# Generalised kinship coefficients
# NB: Only limited patterns supported by Karigl
gKinship_Karigl = function(x, pattern, mem, Xchrom = FALSE, debug = FALSE) {
  if(isDeterministic(pattern))
    stop2("Method 'K' does not support deterministic kinship patterns")

  type = match(paste(lengths(pattern), collapse ="-"),
               c("1", "2", "3", "4", "2-2"))
  if(is.na(type))
    stop2("Method 'K' does not support this pattern")

  if(type == 5 && isDistinct(pattern))
    stop2("The `K` method requires non-distinct blocks. Received: ", gip2string(pattern))

  als = unlist(pattern, use.names = FALSE)

  switch(type,
    1,
    phi2(als[1], als[2], X = Xchrom, mem = mem),
    phi3(als[1], als[2], als[3], X = Xchrom, mem = mem),
    phi4(als[1], als[2], als[3], als[4], X = Xchrom, mem = mem),
    phi22(als[1], als[2], als[3], als[4], X = Xchrom, mem = mem)
  )
}



# Recursion functions -----------------------------------------------------

phi2 = function(a, b, X = FALSE, mem) { # X irrelevant
  mem$i2 = mem$i2 + 1
  if(a*b == 0)
    return(0)
  mem$KIN[[a, b]]
}

phi3 = function(a, b, c, X = FALSE, mem) {
  mem$i3 = mem$i3 + 1
  if(a*b*c == 0)
    return(0)

  ANC = mem$ANC
  if(!ANC[a,b] || !ANC[b,c] || !ANC[a,c])
    return(0)

  # Sort: a >= b >= c
  if(a < b) {tmp = a; a = b; b = tmp}
  if(b < c) {tmp = b; b = c; c = tmp}
  if(a < b) {tmp = a; a = b; b = tmp}

  # Lookup in array
  if(!is.na(res <- mem$KIN3[[a,b,c]]))
    return(res)

  ### Recurse (assumes a,b,c sorted)
  mem$i3r = mem$i3r + 1
  FIDX = mem$FIDX
  MIDX = mem$MIDX

  if(X && mem$SEX[a] == 1) {
    if(a == b && a == c)
      res = 1
    else if(a == b)
      res = phi2(MIDX[a], c, X = X, mem = mem)
    else
      res = phi3(MIDX[a], b, c, X = X, mem = mem)
  }
  else {
    if(a == b && a == c)
      res = (1 + 3*phi2(FIDX[a], MIDX[a], X = X, mem = mem))/4
    else if(a == b)
      res = (phi2(a, c, X = X, mem = mem) +
               phi3(FIDX[a], MIDX[a], c, X = X, mem = mem))/2
    else
      res = (phi3(FIDX[a], b, c, X = X, mem = mem) +
               phi3(MIDX[a], b, c, X = X, mem = mem))/2
  }

  mem$KIN3[[a,b,c]] = res
  res
}

phi4 = function(a, b, c, d, X = FALSE, mem) {
  mem$i4 = mem$i4 + 1
  if(a*b*c*d == 0)
    return(0)

  ANC = mem$ANC
  if(!ANC[a,b] || !ANC[b,c] || !ANC[c,d] || !ANC[b,d] || !ANC[a,c] || !ANC[a,d])
    return(0)

  # Sort: a >= b >= c >= d
  if(a < b) {tmp = a; a = b; b = tmp}
  if(a < c) {tmp = a; a = c; c = tmp}
  if(a < d) {tmp = a; a = d; d = tmp}
  if(b < c) {tmp = b; b = c; c = tmp}
  if(b < d) {tmp = b; b = d; d = tmp}
  if(c < d) {tmp = c; c = d; d = tmp}

  # Lookup in array
  if(!is.na(res <- mem$KIN4[[a,b,c,d]]))
    return(res)

  ### Recurse (assumes a,b,c,d sorted)
  mem$i4r = mem$i4r + 1
  FIDX = mem$FIDX
  MIDX = mem$MIDX

  if(X && mem$SEX[a] == 1) {
    if(a == b && a == c && a == d)
      res = 1
    else if(a == b && a == c)
      res = phi2(a, d, X = X, mem = mem)
    else if(a == b)
      res = phi3(MIDX[a], c, d, X = X, mem = mem)
    else
      res = phi4(MIDX[a], b, c, d, X = X, mem = mem)
  }
  else {
    if(a == b && a == c && a == d)
      res = (1 + 7*phi2(FIDX[a], MIDX[a], X = X, mem = mem))/8
    else if(a == b && a == c)
      res = (phi2(a, d, X = X, mem = mem) +
               3*phi3(FIDX[a], MIDX[a], d, X = X, mem = mem))/4
    else if(a == b)
      res = (phi3(a, c, d, X = X, mem = mem) +
               phi4(FIDX[a], MIDX[a], c, d, X = X, mem = mem))/2
    else
      res = (phi4(FIDX[a], b, c, d, X = X, mem = mem) +
               phi4(MIDX[a], b, c, d, X = X, mem = mem))/2
  }

  mem$KIN4[[a,b,c,d]] = res
  res
}

phi22 = function(a, b, c, d, X = FALSE, mem = NULL) {#cat(a,b,c,d, "\n")
  mem$i22 = mem$i22 + 1
  if(a*b*c*d == 0)
    return(0)

  ANC = mem$ANC
  if(!ANC[a,b] || !ANC[c,d])
    return(0)

  # If no between-pair relations, factorise
  if(!ANC[a,c] && !ANC[a,d] && !ANC[b,c] && !ANC[b,d])
    return(phi2(a, b, X = X, mem = mem) * phi2(c, d, X = X, mem = mem))

  # Sort: a >= b,c,d; c >= d; if(a == c) then b >= d
  s = c(min(a,b), max(a,b), min(c,d), max(c,d)) # d,c,b,a
  if(s[4] < s[2] || (s[4] == s[2] && s[3] < s[1]))
    s[] = s[c(3,4,1,2)]
  a = s[4]; b = s[3]; c = s[2]; d = s[1]

  # Lookup in array
  if(!is.na(res <- mem$KIN22[[a,b,c,d]]))
    return(res)

  ### Recurse (assumes a,b,c,d sorted)
  mem$i22r = mem$i22r + 1
  FIDX = mem$FIDX
  MIDX = mem$MIDX

  if(X && mem$SEX[a] == 1) {
    if(a == b && a == c && a == d)
      res = 1
    else if(a == b && a == c)
      res = phi2(MIDX[a], d, X = X, mem = mem)
    else if(a == b)
      res = phi2(c, d, X = X, mem = mem)
    else if(a == c)
      res = phi3(MIDX[a], b, d, X = X, mem = mem)
    else
      res = phi22(MIDX[a], b, c, d, X = X, mem = mem)
  }
  else {
    if(a == b && a == c && a == d)
      res = (1 + 3*phi2(FIDX[a], MIDX[a], X = X, mem = mem))/4
    else if(a == b && a == c)
      res = (phi2(a, d, X = X, mem = mem) +
               phi3(FIDX[a], MIDX[a], d, X = X, mem = mem))/2
    else if(a == b) { #NB modification to allow inbred founders!
      if(mem$isFounder[a])
        res = 0.5*phi2(c, d, X = X, mem = mem) * (1 + mem$INB[a])
      else
        res = (phi2(c, d, X = X, mem = mem) +
                 phi22(FIDX[a], MIDX[a], c, d, X = X, mem = mem))/2
    }
    else if(a == c)
      res = (2*phi3(a, b, d, X = X, mem = mem) +
               phi22(FIDX[a], b, MIDX[a], d, X = X, mem = mem) +
               phi22(MIDX[a], b, FIDX[a], d, X = X, mem = mem))/4
    else
      res = (phi22(FIDX[a], b, c, d, X = X, mem = mem) +
               phi22(MIDX[a], b, c, d, X = X, mem = mem))/2
  }

  mem$KIN22[[a,b,c,d]] = res
  res
}
