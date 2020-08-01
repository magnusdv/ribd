#' Two-locus kinship coefficients
#'
#' Computes the two-locus kinship coefficient of a pair of pedigree members, at
#' a given recombination rate.
#'
#' Let A, B be two pedigree members, and L1, L2 two loci with a given
#' recombination rate rho. The two-locus kinship coefficient
#' \eqn{\phi_{AB}(rho)} is defined as the probability that random gametes
#' segregating from A and B has IBD alleles at both L1 and L2 simultaneously.
#'
#' The implementation is based on the recursive algorithm described by Thompson
#' (1988).
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object.
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param rho A numeric vector of recombination rates; all entries must be in
#'   the interval \eqn{[0, 0.5]}.
#' @param recombinants A logical of length 2, applicable only when `ids` has
#'   length 2. When given, it indicates whether each of the two gametes is a
#'   recombinant or non-recombinant. This parameter is mainly used by
#'   [twoLocusIBD()].
#' @param verbose A logical.
#' @param debug A logical. If TRUE, detailed messages are printed during the
#'   recursion process.
#'
#' @references E. A. Thompson (1988). _Two-locus and Three-locus Gene Identity
#'   by Descent in Pedigrees_. IMA Journal of Mathematics Applied in Medicine &
#'   Biology, vol. 5.
#'
#'
#' @examples
#' ######################
#' # Example 1: Full sibs
#' ######################
#' x = nuclearPed(2)
#'
#' k_0 = twoLocusKinship(x, ids = 3:4, rho = 0)
#' k_0.5 = twoLocusKinship(x, ids = 3:4, rho = 0.5)
#'
#' stopifnot(k_0 == 1/4, k_0.5 == 1/16)
#'
#'
#' ##################################################
#' # Example 2: Reproducing Fig. 3 in Thompson (1988)
#' # Note that in the article, curve (a) is wrong.
#' # See Erratum: https://doi.org/10.1093/imammb/6.1.1
#' ##################################################
#'
#' # Pedigrees (a) - (d)
#' ped.a = linearPed(3)
#' ped.b = halfCousinPed(0, removal = 1)
#' ped.c = cousinPed(1)
#' ped.d = doubleCousins(1, 1, half1 = TRUE, half2 = TRUE)
#'
#' peds = list(
#'   a = list(ped = ped.a, ids = c(1,7)),
#'   b = list(ped = ped.b, ids = leaves(ped.b)),
#'   c = list(ped = ped.c, ids = leaves(ped.c)),
#'   d = list(ped = ped.d, ids = leaves(ped.d))
#' )
#'
#' # Recombination values
#' rseq = seq(0, 0.5, length = 20)
#'
#' # Compute two-locus kinship coefficients
#' kvals = sapply(peds, function(x) twoLocusKinship(x$ped, x$ids, rseq))
#'
#' # Plot
#' matplot(rseq, kvals, type = "l", lwd = 2)
#' legend("topright", names(peds), col = 1:4, lty = 1:4, lwd = 2)
#'
#' @importFrom utils combn
#' @export
twoLocusKinship = function(x, ids, rho, recombinants = NULL, verbose = FALSE, debug = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) < 2) stop2("Argument `ids` must have length at least 2")
  if(!is.numeric(rho))
    stop2("Argument `rho` must be numeric")
  if(any(rho < 0 | rho > 0.5))
    stop2("Argument `rho` cannot have entries outside the interval [0, 0.5]: ", rho[rho < 0 | rho > 0.5])
  if(!is.null(recombinants) && length(ids) > 2)
    stop2("Argument `recombinants` must be NULL when `ids` has length > 2")
  if(debug && (length(ids) > 2 || length(rho) > 1))
    stop2("Debugging mode is only allowed when `ids` has length 2, and `rho` has length 1")

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  ids_int = internalID(x, ids)

  # Convert recombination conditions from logical to list(r = , nr = )
  if(!is.null(recombinants))
    rList = list(nr = ids[!recombinants], r = ids[recombinants])
  else
    rList = NULL

  # Simplest case: Allows verbose output and/or debugging mode
  if(length(ids) == 2 && length(rho) == 1) {
    A = C = c(ids_int[1], -1) # using negative numbers to enforce independent gametes
    B = D = c(ids_int[2], -2)

    # Counters for verbose output
    counters = c("i","ilook","irec","eq7","eq8","eq9a","eq9b","eq10","eq11a","eq11b")

    # Initialise memory storage
    mem = initialiseTwoLocusMemo(x, rho = rho, recomb = rList, counters = counters)

    # Compute!
    phi11 = twoLocKin(A, B, C, D, mem, indent = ifelse(debug, 0, NA))

    # Print info
    if(verbose) {
      init_time = format(mem$initTime, digits = 4)
      tot_time = format(Sys.time()-mem$st, digits = 4)
      message(glue::glue("
        Calls = {mem$i}
        Lookups = {mem$ilook}
        Recursions: {mem$irec}
          eq. 7  : {mem$eq7}
          eq. 8  : {mem$eq8}
          eq. 9a : {mem$eq9a}
          eq. 9b : {mem$eq9b}
          eq. 10 : {mem$eq10}
          eq. 11a: {mem$eq11a}
          eq. 11b: {mem$eq11b}
        Total time used: {tot_time}"))
    }

    return(phi11)
  }

  # Storage template: Reset for each rho
  memTemplate = as.list(initialiseTwoLocusMemo(x, rho = NULL, recomb = rList))

  # If single pair: Return only vector of coefficients
  if(length(ids) == 2) {
    A = C = c(ids_int[1], -1) # using negative numbers to enforce independent gametes
    B = D = c(ids_int[2], -2)

    coefs = vapply(rho, function(r) {
      mem = as.environment(memTemplate)
      mem$rho = r
      twoLocKin(A, B, C, D, mem, indent = NA)
    }, FUN.VALUE = 0)

    return(coefs)
  }

  # If length(ids) > 2: Do all unordered pairs; return data.frame

  pairs = combn(ids_int, 2, simplify = FALSE)
  pairs = c(pairs, lapply(seq_along(ids_int), function(i) c(i,i)))

  coefs = lapply(rho, function(r) {
    mem = as.environment(memTemplate)
    mem$rho = r

    unlist(lapply(pairs, function(p) {
      A = C = c(p[1], -1)
      B = D = c(p[2], -2)
      twoLocKin(A, B, C, D, mem, indent = NA)
    }))
  })

  # Build result data frame
  labs = labels(x)
  idcols = do.call(rbind, pairs)
  idcols[] = labs[idcols]
  res = data.frame(id1 = idcols[,1],
                   id2 = idcols[,2],
                   rho = rep(rho, each = length(pairs)),
                   phi2 = unlist(coefs),
                   stringsAsFactors = FALSE)

  res
}

# Thompson (1988) k2(J(A,B), L(C,D)). A = c(id, segind) a.s.o.
twoLocusK2 = function(x, J, L, rho, verbose = FALSE, debug = FALSE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")
  if(!is.list(J) || length(J) != 2)
    stop2("Argument `J` must be a list of length 2: ", J)
  if(!is.list(L) || length(L) != 2)
    stop2("Argument `L` must be a list of length 2: ", L)
  if(!is.numeric(rho) && length(rho) != 1)
    stop2("Argument `rho` must be a single numeric")

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  # Counters for verbose output
  counters = c("i","ilook","irec","eq7","eq8","eq9a","eq9b","eq10","eq11a","eq11b")

  mem = initialiseTwoLocusMemo(x, rho = rho, counters = counters)

  # Convert to internal numeric labels
  convertAndFix = function(S, segind) {
    S[1] = internalID(x, S[1])
    if(length(S) == 1)  # If segregation index missing, insert supplied value
      S = c(S, segind)
    else if(S[2] > 0)
      S[2] = internalID(x, S[2])
    S
  }

  # Taken (negative) segreg indices
  takenSeg = min(0, J[[1]][2], J[[2]][2], L[[1]][2], L[[2]][2])

  A = convertAndFix(J[[1]], takenSeg - 1)
  B = convertAndFix(J[[2]], takenSeg - 2)
  C = convertAndFix(L[[1]], takenSeg - 3)
  D = convertAndFix(L[[2]], takenSeg - 4)

  twoLocKin(A, B, C, D, mem, indent = ifelse(debug, 0, NA))
}




####################
# Internal functions
####################

twoLocKin = function(A, B, C, D, mem, indent = 0) {
  mem$i = mem$i + 1
  # Each A, B, C, D is a parent-child pair

  # If any founders -> 0
  if(A[1]*B[1]*C[1]*D[1] == 0) {message("Outside of pedigree! Can this happen?"); return(0)}

  # Sort: a >= b,c,d; c >= d; if(a == c) then b >= d
  plist = sortPairs(A,B,C,D)
  printMess(plist, indent)

  A = plist[[1]]; B = plist[[2]]; C = plist[[3]]; D = plist[[4]]
  a = A[1]; b = B[1]; c = C[1]; d = D[1]

  k2_recurse = function(A,B,C,D) {
    mem$irec = mem$irec + 1
    a = A[1]; b = B[1]; c = C[1]; d = D[1]

    if(identical(A, B))
      stop2("Identical gametes at first locus; not implemented. (This should not occur in a two-locus kinship computation.)")

    if(a > b && a > c) { # eq. 7
      recurse_eq7(A,B,C,D,mem,indent)
    }
    else if(a == b && a > c) { # eq. 8
      recurse_eq8(A,B,C,D,mem,indent)
    }
    else if(a > b && a == c) { #eq. 9a or 9b
      if(A[2] == C[2])
        recurse_eq9a(A,B,C,D,mem,indent)
      else
        recurse_eq9b(A,B,C,D,mem,indent)
    }
    else if(a == b && a == c && a > d) { # eq. 10
      recurse_eq10(A,B,C,D,mem,indent)
    }
    else if((a == b && a == c && a == d)) {
      # This case needs further branching!

      if(A[2] == C[2] && B[2] == D[2]) # k2(J(A1, A2), L(A1, A2))
        recurse_eq11a(A,B,C,D,mem,indent)

      else  # k2(J(A1, A2), L(A1, A3)) or k2(J(A1, A2), L(A3, A4))
        recurse_eq11b(A,B,C,D,mem,indent)
    }
  }

  isFou = mem$isFounder[a]
  ANC = mem$anc

  # If no common ancestors, return 0
  if(!(ANC[a,b] && ANC[c,d]))
    return(printAndReturn(0, indent))

  # TODO: Reduce dimensionality if possible, as in weeks & lange
  # TODO: Use symmetry if (A1, B2, A1, B2), as in Thompson
  # TODO: Special cases not covered in conditional recomb/nonrecomb: Indicate segreg!

  # Lookup in array; compute if necessary.
  mem$ilook = mem$ilook + 1
  res = mem$k2[[toString(plist)]]
  if(is.null(res))
    res = mem$k2[[toString(plist)]] = k2_recurse(A, B, C, D)

  printAndReturn(res, indent)
}


printMess = function(plist, indent) {
  if(is.na(indent)) return()
  pp = sapply(plist, function(p) paste(p, collapse = ":"))
  message(sprintf("%sJ(%s, %s), L(%s, %s)",
                  strrep(" ", indent), pp[1], pp[2], pp[3], pp[4]))
}

printAndReturn = function(res, indent, comment = NULL) {
  if(!is.na(indent))
    message(strrep(" ", indent), res, comment)
  res
}


### The recursions (eqs. 7-11 in Thompson (1988) + additional cases in Weeks&Lange)

recurse_eq7 = function(A,B,C,D,mem,indent) {
  mem$eq7 = mem$eq7 + 1
  a = A[1]
  if(mem$isFounder[a])
    return(0)

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]

  t1 = twoLocKin(c(FF, a), B, C, D, mem, indent = indent + 2)
  t2 = twoLocKin(c(MM, a), B, C, D, mem, indent = indent + 2)
  0.5 * (t1 + t2)
}

recurse_eq8 = function(A,B,C,D,mem,indent) {
  mem$eq8 = mem$eq8 + 1
  a = A[1]
  phi_cd = mem$k1[[C[1], D[1]]] # kinship of c and d

  if(mem$isFounder[a]) {
    if(mem$isCompletelyInbred[a])
      return(phi_cd)
    else
      return(0.5 *  phi_cd)
  }

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]
  0.5 * (phi_cd + twoLocKin(c(FF, a), c(MM, a+100), C, D, mem, indent = indent + 2))
}

recurse_eq9a = function(A,B,C,D,mem,indent) {
  mem$eq9a = mem$eq9a + 1
  a = A[1]
  if(mem$isFounder[a])
    return(0)

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]
  rho = mem$rho


  # Terms used in various formulas
  # NB (a+100): used to separate meioses involved in selfing, i.e. if FF = MM.
  t1 = function() twoLocKin(c(FF, a), B, c(FF, a), D, mem, indent = indent + 2)
  t2 = function() twoLocKin(c(MM, a), B, c(MM, a), D, mem, indent = indent + 2)
  t3 = function() twoLocKin(c(FF, a), B, c(MM, a+100), D, mem, indent = indent + 2)
  t4 = function() twoLocKin(c(MM, a), B, c(FF, a+100), D, mem, indent = indent + 2)

  if(mem$nonrecomb[a])              # force non-rec
    res = 0.5 * (1-rho) * (t1() + t2())
  else if(mem$recomb[a])          # force rec
    res = 0.5 * rho * (t3() + t4())
  else                                  # normal recursion
    res = 0.5 * ((1-rho) * (t1() + t2()) + rho * (t3() + t4()))

  res
}

recurse_eq9b = function(A,B,C,D,mem,indent) {
  mem$eq9b = mem$eq9b + 1
  a = A[1]

  if(mem$isFounder[a])
    return(0)

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]

  # NB (a+100): used to separate meioses involved in selfing, i.e. if FF = MM.
  t1 = twoLocKin(c(FF, a), B, c(FF, a), D, mem, indent = indent + 2)
  t2 = twoLocKin(c(MM, a), B, c(MM, a), D, mem, indent = indent + 2)
  t3 = twoLocKin(c(FF, a), B, c(MM, a+100), D, mem, indent = indent + 2)
  t4 = twoLocKin(c(MM, a), B, c(FF, a+100), D, mem, indent = indent + 2)

  0.25 * (t1 + t2 + t3 + t4)
}

recurse_eq10 = function(A,B,C,D,mem,indent) {
  mem$eq10 = mem$eq10 + 1
  a = A[1]

  if(mem$isFounder[a])
    return(0)

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]
  k1 = mem$k1
  rho = mem$rho
  d = D[1]

  t3 = twoLocKin(c(FF, a), c(MM, a), c(FF, a+100), D, mem, indent = indent + 2)
  t4 = twoLocKin(c(FF, a), c(MM, a), c(MM, a+100), D, mem, indent = indent + 2)

  if(A[2] == C[2]) { # EAT eq 10
    s = 0.25 * (k1[[FF, d]] + k1[[MM, d]] + t3 + t4)

    if(mem$nonrecomb[a])
      res = (1-rho) * s
    else if(mem$recomb[a])
      res = rho * s
    else
      res = s
  }
  else {
    if(mem$nonrecomb[a] || mem$recomb[a])
      stop2("eq 10b conditional: Not implemented")
    res = 0.25 * (2*k1[[a, d]] + t3 + t4)
  }

  res
}


recurse_eq11a = function(A,B,C,D,mem,indent) { # Case k2(A1,A2; A1,A2): Eq. 11 in EAT
  mem$eq11a = mem$eq11a + 1
  a = A[1]
  forceNonRec = mem$nonrecomb[a]
  forceRec = mem$recomb[a]
  rho = mem$rho
  R = .5 * ((1-rho)^2 + rho^2)

  if(mem$isFounder[a]) {
    if(mem$isCompletelyInbred[a])
      return(1)

    if(forceNonRec && forceRec)
      res = 0
    else if(forceNonRec)
      res = .5 * (1-rho)^2
    else if(forceRec)
      res = .5 * rho^2
    else
      res = R

    return(res)
  }

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]
  k1 = mem$k1

  if(forceNonRec && forceRec) {
    res = rho*(1-rho)*k1[[MM, FF]] # without the factor two (either R-NR or NR-R)
  } else {
    t4 = twoLocKin(c(MM, a), c(FF, a+100), c(MM, a), c(FF, a+100), mem, indent = indent + 2)
    if(forceNonRec)
      res = .5 * (1-rho)^2 * (1 + t4)
    else if(forceRec)
      res = .5 * rho^2 * (1 + t4)
    else
      res = 2*rho*(1-rho)*k1[[MM, FF]] + R*(1 + t4)
    }

  res
}


recurse_eq11b = function(A,B,C,D,mem,indent) { # k2(A1,A2; A1,A3) or k2(A1,A2; A3,A4)
  # Recursion formula given by Weeks&Lange
  mem$eq11b = mem$eq11b + 1
  a = A[1]
  rho = mem$rho

  # Condition on recomb status: only for k2(A1,A2; A1,A3)
  nonrecomb = mem$nonrecomb[a] && A[2] == C[2] && A[2] %in% c(-1, -2)
  recomb = mem$recomb[a] && A[2] == C[2] && A[2] %in% c(-1, -2)
  if(nonrecomb && recomb)
    stop2("Recursion 11b) Special case not implemented")


  if(mem$isFounder[a]) {
    if(nonrecomb) return((1 - rho)/4)
    if(recomb) return(rho/4)
    return(1/4)
  }

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]
  k1 = mem$k1

  parents2L = twoLocKin(c(FF, a), c(MM, a+100), c(FF, a), c(MM, a+100), mem, indent = indent + 2)
  res = 1/4 + 1/2 * k1[[FF, MM]] + 1/4 * parents2L

  if(nonrecomb) res = (1 - rho) * res
  if(recomb) res = rho * res

  res
}


sortPairs = function(A,B,C,D) {
  plist = list(A,B,C,D)

  ### Sort: a >= b,c,d; c >= d; if(a == c) then b >= d
  if (plist[[1]][1] < plist[[2]][1])
    plist[1:2] = plist[2:1]
  if (plist[[3]][1] < plist[[4]][1])
    plist[3:4] = plist[4:3]
  if (plist[[1]][1] < plist[[3]][1] || (plist[[1]][1] == plist[[3]][1] && plist[[2]][1] < plist[[4]][1]))
    plist[] = plist[c(3,4,1,2)]

  ### Ensure segregation indicators are sorted appropriately
  A = plist[[1]]; B = plist[[2]]; C = plist[[3]]; D = plist[[4]]

  # Only necessary when a = b = c
  if (A[1] == B[1] && A[1] == C[1]) {
    if (A[1] == D[1]) { # case 1: All equal
      if(A[2] != C[2] && A[2] == D[2])
        plist[3:4] = plist[4:3]
      else if(A[2] != C[2] && B[2] == C[2])
        plist[1:2] = plist[2:1]
      else if(A[2] != C[2] && B[2] == D[2])
        plist[] = plist[c(2:1,4:3)]
    }
    else { # case 2: a = b = c != d
      if(A[2] != C[2] && B[2] == C[2])
        plist[1:2] = plist[2:1]
    }
  }

  plist
}


