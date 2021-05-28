#' Multi-person IBD coefficients
#'
#' Computes the probabilities (coefficients) of all possible patterns of
#' identity by descent (IBD) sharing at a single locus, among N>1 non-inbred
#' members of a pedigree. The reported coefficients are "condensed" in the sense
#' that allele ordering within each individual is ignored. For N = 2, the result
#' should agree with the traditional "kappa" coefficients, as computed by
#' [kappaIBD()]. This function is under development, and should be regarded as
#' experimental. For now, the only cases handled are those with: N = 2 or 3,
#' autosomal locus.
#'
#' Consider N members of a pedigree, i1, i2, ... iN.  A pattern of IBD sharing
#' between these individuals is a sequence of N ordered pairs of labels, (a1_1,
#' a1_2), (a2_1, a2_2), ... (aN_1, aN_2), where ai_1 and ai_2 represent the
#' paternal and maternal allele of individual i, respectively. Equality of
#' labels means that the corresponding alleles are IBD, and vice versa.
#'
#' We say that two IBD patterns are equivalent if one can be transformed into
#' the other by some combination of
#'
#' * renaming the labels (without changing the structure)
#'
#' * swapping the paternal/maternal labels of some individuals
#'
#' Each equivalence class has a "minimal" element, using integer labels, and
#' being minimal with respect to standard sorting. For example, the minimal
#' element equivalent to `(a,c),(d,c),(b,b)` is `(1,2),(2,3),(4,4)`.
#'
#' @param x A `ped` object.
#' @param ids A vector of ID labels.
#' @param complete A logical. If FALSE, only IBD patterns with nonzero
#'   probability are included in the output.
#' @param verbose A logical. If TRUE, some computational details are printed.
#'
#' @return A data frame in which each row corresponds to an equivalence class of
#'   multi-person IBD patterns. The first column gives the calculated
#'   probability, followed by one column for each `ids` individual, describing
#'   the minimal element of the equivalence class. (See Details.) If `complete =
#'   FALSE` (the default) rows with probability 0 are removed.
#'
#' @examples
#' ### Trivial example: Trio ###
#' x = nuclearPed(1)
#' ids = 1:3
#' multiPersonIBD(x, ids, complete = TRUE)
#'
#' ### Example due to Peter Green ###
#' # Three (pariwise) cousins arranged in two different ways,
#' # with different 3-way IBD coefficients.
#'
#' threeCousins1 = ped(
#'   id  = c('gf','gm','gf1','gf2','gf3','gm1','gm2','gm3',
#'           'f1','f2','f3','m1','m2','m3','c1','c2','c3'),
#'   fid = c(0,0,0,0,0,0,0,0,'gf1','gf2','gf3','gf','gf','gf',
#'           'f1','f2','f3'),
#'   mid = c(0,0,0,0,0,0,0,0,'gm1','gm2','gm3','gm','gm','gm',
#'           'm1','m2','m3'),
#'   sex = c(1,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1))
#'
#' threeCousins2 = ped(
#'   id  = c('gf1','gf2','gf3','gm1','gm2','gm3','f1','f2','f3',
#'           'm1','m2','m3','c1','c2','c3'),
#'   fid = c(0,0,0,0,0,0,'gf2','gf3','gf1','gf3','gf1','gf2',
#'           'f1','f2','f3'),
#'   mid = c(0,0,0,0,0,0,'gm2','gm3','gm1','gm3','gm1','gm2',
#'           'm1','m2','m3'),
#'   sex = c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1))
#'
#' ids = c('c1','c2','c3')
#' multiPersonIBD(threeCousins1, ids)
#' multiPersonIBD(threeCousins2, ids)
#'
#' @export
multiPersonIBD = function(x, ids, complete = FALSE, verbose = FALSE) {
  N = length(ids)
  if(N < 2)
    stop2("`length(ids)` must be at least 2")
  if(N > 6)
    stop2("length(ids) > 6: Not implemented yet")

  # Stop if any of `ids` are inbred
  inb = inbreeding(x, ids)
  if(any(isInbred <- inb > .Machine$double.eps)) {
    message(paste0(sprintf(" Individual '%s' is inbred (f = %g)",
                           ids[isInbred], inb[isInbred]), collapse = "\n"))
    stop2("This function requires all `ids` individuals to be non-inbred")
  }

  # Pedigree order: Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  x = foundersFirst(x)

  # Setup memoisation
  mem = initialiseGKMemo(x, counters = c("i", "itriv", "iimp", "ifound", "ilook", "irec"))

  allPatterns = MULTIPATTERNS_NONINBRED[[N]]
  usePatterns = removeImpossiblePatterns(allPatterns, x, ids, verbose = verbose)

  # Vector of `ids` repeated twice (needed in the generalised below)
  ids2 = rep(ids, each = 2)

  coefs = apply(usePatterns, 1, function(r) {
    kp = r[1:(2*N)]
    weight = r[2*N + 1]
    g = generalisedKinship(x, split(ids2, kp), mem = mem)
    2^N * weight * g
  })

  # Print computational summary
  if(verbose)
    printCounts2(mem)

  # Collect into data frame
  glist = lapply(seq(1, 2*N, by = 2), function(i)
    paste(usePatterns[, i], usePatterns[, i+1], sep = " "))

  res = as.data.frame(glist, col.names = as.character(ids), optional = TRUE,
                      stringsAsFactors = FALSE)
  res = cbind(Prob = coefs, res)

  # Keep only nonzero rows
  if(!complete)
    res = res[res$Prob > 0, , drop = FALSE]

  # Return
  res
}


#### Utility functions used in multiPersonIBD ####

# Reduce a sequence of alleles to "standard" form
# a,c,a,d     -> 1,2,1,3
# 2,1,4,3,1,1 -> 1,2,3,4,2,2
# NB: simpler but slower: as.integer(factor(x, levels = unique.default(x)))
standardAlleleSeq = function(x) { # x a vector of even length
  res = integer(length(x))
  a = x[1]
  i = 1L
  while(TRUE) {
    res[x == a] = i
    a = x[match(0, res)]
    if(is.na(a))
      break
    i = i + 1
  }
  res
}

#' Minimal IBD pattern
#'
#' Compute the minimal form of given multiperson IBD pattern.
#'
#' @param x An integer vector of even length.
#'
#' @return An integer vector of the same length as `x`.
#'
#' @examples
#' v = c(1,2,2,3)
#' stopifnot(identical(minimalPattern(v), c(1,2,1,3)))
#'
#' @export
minimalPattern = function(x) {
  if(is.matrix(x)) {
    res = apply(x, 1, minimalPattern)
    return(res)
  }

  if(!is.numeric(x) || any(as.integer(x) != x))
    stop("The input must be a vector of integers:", toString(x))

  n = length(x)
  if(n %% 2 != 0)
    stop("Input vector must have even length, not ", toString(x))

  # Utility for inserting a certain value at given positions.
  # If "value" already exist in the vector, these are replaced with a larger dummy
  insertAt = function(y, pos, value) {
    y[y == value] = max(y) + 1 # dummy
    y[pos] = value
    y
  }

  nPairs = n/2
  even = 2 * seq_len(nPairs)
  odd = even - 1
  i = 0

  # Loop through pairs of alleles (i.e. genotypes)
  for(k in seq_len(nPairs)) {
    a = x[2*k - 1]
    b = x[2*k]

    if(a == b) {
      if(a > i) {
        x = insertAt(x, x == a, i + 1)
        i = i + 1
      }
      next
    }

    # How many alleles are not seen before?
    n.new = sum(a > i, b > i)

    if(n.new == 2) { # both new!
      idx.a = x == a
      idx.b = x == b
      grps.a = idx.a[odd] | idx.a[even]
      grps.b = idx.b[odd] | idx.b[even]

      # Swap labels if the first genotype with exactly one of a,b has b.
      diffs = grps.a != grps.b
      swap = any(diffs) && grps.b[match(TRUE, diffs)]

      x = insertAt(x, idx.a, value = if(swap) i+2 else i+1)
      x = insertAt(x, idx.b, value = if(swap) i+1 else i+2)
      i = i + 2
    }
    else if (n.new == 1) { # only one new
      x = insertAt(x, x == max(a,b), i + 1)
      i = i + 1
    }

    # Always ensure the resulting pair is sorted
    if((g1 <- x[2*k - 1]) > (g2 <- x[2*k])) {
      x[2*k - 1] = g2
      x[2*k] = g1
    }
  }

  # Return modified vector
  x
}

# Produce all permutations made by swapping the two alleles within (not between) individuals
expandSwaps = function(x, makeUnique = TRUE) { # x a vector of even length
  n = length(x)/2

  res = matrix(0L, nrow = 2*n, ncol = 2^n)
  for(i in 1:n) {
    rows = c(2*i-1, 2*i)
    g = x[rows]
    swap = rep(c(FALSE, TRUE), each = 2^(n-i), length.out = 2^n)
    res[rows, !swap] = g
    res[rows, swap] = g[2:1]
  }

  res = apply(res, 2, standardAlleleSeq)

  if(makeUnique)
    res = unique.matrix(res, MARGIN = 2)

  t.default(res)
}


# Remove impossible patterns based on pairwise kappa
removeImpossiblePatterns = function(patterns, x, ids, verbose = TRUE) {
  N = length(ids)
  if(length(ids) < 3) return(patterns[, 1:(2*N + 1)])

  # Original number of patterns
  nr = nrow(patterns)

  kappas = kappaIBD(x, ids)

  nonz = lapply(1:nrow(kappas), function(i) (0:2)[kappas[i, 3:5] > 0])
  names(nonz) = sapply(combn(N, 2, simplify = FALSE), paste, collapse = "-")
  for(p in names(nonz)) {
    goodrows = patterns[, p] %in% nonz[[p]]
    patterns = patterns[goodrows, , drop = FALSE]
  }

  # Report change
  if(verbose)
    message(sprintf("State space reduction: %d --> %d", nr, nrow(patterns)))

  patterns
}
