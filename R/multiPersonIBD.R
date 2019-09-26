#' Multi-person IBD coefficients
#'
#' Computes the probabilities (coefficients) of all possible patterns of
#' identity by descent (IBD) sharing at a single locus, among N>1 non-inbred
#' members of a pedigree. For now, only N = 2 and 3, and only autosomal loci are
#' implemented. The reported coefficients are "condensed" in the sense that
#' allele ordering within each individual is ignored. For N = 2, the result
#' should agree with the traditional "kappa" coefficients, as computed by
#' [kappaIBD()].
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
#' @export
#'
#' @examples
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
multiPersonIBD = function(x, ids, complete = F, verbose = F) {
  N = length(ids)
  if(N < 2)
    stop2("`length(ids)` must be at least 2")
  if(N > 3)
    stop2("length(ids) > 3: Not implemented yet")

  # Stop if any of `ids` are inbred
  inb = inbreeding(x)[ids]
  if(any(isInbred <- inb > .Machine$double.eps))
    stop2(paste0(c(sprintf(" Individual '%s' is inbred (f = %g)",
                           ids[isInbred], inb[isInbred]),
                   "This function requires non-inbred individuals."), collapse = "\n"))

  # If founders: add parents
  for(i in intersect(ids, founders(x)))
    x = addParents(x, i, verbose = F)

  # Pedigree order: force founders first (TODO: is this needed here?)
  x = foundersFirst(x)

  # Pedigree order: Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  # Setup memoisation
  mem = initialiseGKMemo(x, counters = c("i", "itriv", "iimp", "ifound", "ilook", "irec"))

  #
  if(N == 2) {
    allPatterns = matrix(c(
      1,2,1,2, #(2)
      1,2,1,3, #(1)
      1,2,3,4),#(0)
      byrow = T, ncol = 4)
  }
  else if(N == 3) {
    # All 16 patterns of 3 indivs, in standard form
    allPatterns = matrix(c(
      1,2,1,2,1,2, #(2,2,2)
      1,2,1,2,1,3, #(2,1,1)
      1,2,1,3,1,2, #(1,2,1)
      1,2,1,3,1,3, #(1,1,2)
      1,2,1,3,1,4, #(1,1,1)
      1,2,1,3,2,3, #(1,1,1)
      1,2,1,2,3,4, #(2,0,0)
      1,2,3,4,1,2, #(0,2,0)
      1,2,3,4,3,4, #(0,0,2)
      1,2,1,3,2,4, #(1,1,0)
      1,2,1,3,3,4, #(1,0,1)
      1,2,3,4,1,3, #(0,1,1)
      1,2,1,3,4,5, #(1,0,0)
      1,2,3,4,1,5, #(0,1,0)
      1,2,3,4,3,5, #(0,0,1)
      1,2,3,4,5,6),#(0,0,0)
      byrow = T, ncol = 6)
  }

  # Vector of father,mother for the `ids` individuals
  famo = unlist(lapply(ids, parents, x=x))

  coefs = apply(allPatterns, 1, function(r) {
    kps = expandSwaps(r, makeUnique = T)
    sum(apply(kps, 1, function(kp)
      generalisedKinship(x, split(famo, kp), mem = mem)))
    })

  # Print computational summary
  if(verbose)
    printCounts2(mem)

  # Collect into data frame
  glist = lapply(seq(1, ncol(allPatterns), by = 2), function(i)
    paste(allPatterns[, i], allPatterns[, i+1], sep=" "))

  res = as.data.frame(glist, col.names = as.character(ids), stringsAsFactors = F)
  res = cbind(Prob = coefs, res)

  # Keep only nonzero rows
  if(!complete)
    res = res[res$Prob > 0, , drop = F]

  # Remove row names
  row.names(res) = NULL

  # Return
  res
}


#### Utility functions used in multiPersonIBD ####

# Reduce a sequence of alleles to "standard" form
# a,c,a,d     -> 1,2,1,3
# 2,1,4,3,1,1 -> 1,2,3,4,2,2
standardAlleleSeq = function(x) { # x a vector of even length
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
expandSwaps = function(x, makeUnique = T) { # x a vector of even length
  n = length(x)/2

  res = matrix(0L, nrow = 2*n, ncol = 2^n)
  for(i in 1:n) {
    rows = c(2*i-1, 2*i)
    g = x[rows]
    swap = rep(c(F,T), each = 2^(n-i), length.out = 2^n)
    res[rows, !swap] = g
    res[rows, swap] = g[2:1]
  }

  res = apply(res, 2, standardAlleleSeq)

  if(makeUnique)
    res = unique.matrix(res, MARGIN = 2)

  t.default(res)
}

