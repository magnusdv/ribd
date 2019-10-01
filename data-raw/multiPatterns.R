# Reduce a sequence of alleles to "standard" form
# a,c,a,d     -> 1,2,1,3
# 2,1,4,3,1,1 -> 1,2,3,4,2,2
# x = a multi-pattern with any labels (a vector of even length)
standardAlleleSeq = function(x) {
  ## slower alternative: as.integer(factor(x, levels = unique.default(x)))
  res = integer(length(x))
  a = x[1]
  i = 1L
  while(T) {
    res[x == a] = i
    a = x[match(0, res)]
    if(is.na(a))
      break
    i = i + 1
  }
  res
}

# Produce all permutations made by swapping the two alleles within (not between) individuals
# x = a multi-pattern (a positive integer vector of even length)
expandSwaps = function(x, makeUnique = T) {
  n = length(x)/2

  res = matrix(0L, nrow = 2*n, ncol = 2^n)
  for(i in 1:n) {
    rows = c(2*i-1, 2*i)
    g = x[rows]
    swap = rep(c(F,T), each = 2^(n-i), length.out = 2^n)
    res[rows, !swap] = g
    res[rows, swap] = g[2:1]
  }

  res = apply(res, 2, function(cc) standardAlleleSeq(cc))

  if(makeUnique)
    res = unique.matrix(res, MARGIN = 2)

  t.default(res)
}

# Find the minimal pattern equivalent to the input.
# x = a multi-pattern (a positive integer vector of even length)
minimalPattern = function(x) {
  n = length(x)
  swaps = expandSwaps(x, makeUnique = F)
  s = apply(swaps, 1, function(v) sum(n^((n-1):0) * (v-1)))
  swaps[which.min(s), ]
}

# Add extra columns with info about each multi-pattern:
# * Number of equivalent "swaps"
# * Inbreeding status for each indiv
# * IBD status for each pair
addExtraInfo = function(x) { # x = matrix where each row is a multi-pattern
  n = ncol(x)/2

  uniqueSwaps = apply(x, 1, function(v) nrow(expandSwaps(v, makeUnique = T)))
  x = cbind(x, swaps = uniqueSwaps)

  # Inbreeding info
  for(i in 1:n) {
    x = cbind(x, x[,2*i-1] == x[,2*i])
  }
  colnames(x)[2*n + 1 + 1:n] = as.character(1:n)

  if(n < 2)
    return(x)

  # IBD status for each pair
  prs = combn(n, 2, simplify = F)
  for(p in prs) {
    i = p[1]; j = p[2]
    i1 = x[,2*i-1]
    i2 = x[,2*i]
    j1 = x[,2*j-1]
    j2 = x[,2*j]
    ibd = as.integer(i1 == j1) + as.integer(i2 == j2) + as.integer(i1 == j2) + as.integer(i2 == j1)
    x = cbind(x, ibd)
  }
  colnames(x)[3*n + 1 + 1:length(prs)] = sapply(prs, paste, collapse = "-")

  x
}
########################

generateNext = function(x, inbred = F) {
  if(is.matrix(x))
    x = lapply(1:nrow(x), function(i) x[i, ])

  res = list()
  for(v in x) {
    mx = max(v)
    for(i in 1:(mx+1)) for(j in (i+1-as.integer(inbred)):max(mx + 1, i+ 1)) {
      w = c(v, i, j)
      wmin = minimalPattern(w)
      res = c(res, list(wmin))
    }
  }
  mat = do.call(rbind, res)
  unique.matrix(mat)
}

L1 = rbind(1:2)
L2 = generateNext(L1)
L3 = generateNext(L2)
L4 = generateNext(L3)
L5 = generateNext(L4)
L6 = generateNext(L5)

MULTIPATTERNS_NONINBRED = lapply(list(L1, L2, L3, L4, L5, L6), addExtraInfo)

usethis::use_data(MULTIPATTERNS_NONINBRED, overwrite = TRUE)
