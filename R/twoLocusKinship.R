#' Two-locus kinship coefficients
#'
#' Computes the two-locus kinship coefficient of a pair of pedigree members, at
#' a given recombination rate.
#'
#' Let A, B be two pedigree members, and L1, L2 two loci with a given
#' recombination rate r. The two-locus kinship coefficient \eqn{\phi_{AB}(r)} is
#' defined as the probability that random gametes segregating from A and B has
#' IBD alleles at both L1 and L2 simultaneously.
#'
#' The implementation is based on the recursive algorithm described by Thompson
#' (1988).
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object.
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param r A number in the interval \eqn{[0, 0.5]}; the recombination rate
#'   between the two loci.
#' @param verbose A logical.
#' @param debug A logical. If TRUE, detailed messages are printed during the
#'   recursion process.
#'
#' @references E. A. Thompson (1988). _Two-locus and Three-locus Gene Identity
#'   by Descent in Pedigrees_ IMA Journal of Mathematics Applied in Medicine &
#'   Biology, vol. 5.
#'
#'
#' @examples
#' ######################
#' # Example 1: Full sibs
#' ######################
#' x = nuclearPed(2)
#'
#' k_0 = twoLocusKinship(x, ids = 3:4, r = 0)
#' k_0.5 = twoLocusKinship(x, ids = 3:4, r = 0.5)
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
#' peds = list(
#'   a = list(ped = linearPed(3), ids = c(1,7)),
#'   b = list(ped = halfCousinPed(0, 1), ids = c(3,7)),
#'   c = list(ped = cousinPed(1), ids = c(5,8)),
#'   d = list(ped = doubleCousins(1, 1, half1 = TRUE, half2 = TRUE), ids = c(5,9))
#' )
#'
#' # Recombination values
#' rseq = seq(0, 0.5, length = 20)
#'
#' # Compute two-locus kinship coefficients
#' kvals = sapply(peds, function(x) twoLocusKinship(x$ped, x$ids, rseq))
#'
#' # Plot
#' matplot(rseq, kvals, type="l", lwd=2, )
#' legend("topright", names(peds), col = 1:4, lty = 1:4)
#'
#' @importFrom utils combn
#' @export
twoLocusKinship = function(x, ids, r, recombinants = NULL, verbose = FALSE, debug = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  # Convert recombination conditions from logical to list(r = , nr = )
  recombList = if(is.null(recombinants)) NULL else list(nr = ids[!recombinants], r = ids[recombinants])

  if(length(ids) == 2) {
    mem = initialiseTwoLocusMemo(x, r = r, recomb = recombList)

    A = C = c(ids_int[1], 0)
    B = D = c(ids_int[2], 0)
    phi11 = twoLocKin(A, B, C, D, mem, indent = ifelse(debug, 0, NA))

    # Print info
    if(verbose) {
      initsecs = sprintf("%.2f", mem$initTime)
      totsecs = sprintf("%.1f", Sys.time()-mem$st)
      print(glue::glue("
                       Calls = {mem$i}
                       Lookups = {mem$ilook}
                       Recursions: {mem$irec}
                          eq. 7: {mem$eq7}
                          eq. 8: {mem$eq8}
                          eq. 9a: {mem$eq9a}
                          eq. 9b: {mem$eq9b}
                          eq. 10: {mem$eq10}
                          eq. 11: {mem$eq11}
                       Total time used: {totsecs} seconds"))
    }

    return(phi11)
  }

  # If length(ids) > 2: Do all unordered pairs; return data.frame
  mem = NULL
  memTemplate = as.list(initialiseTwoLocusMemo(x, r = NULL, recomb = recombList))

  pairs = combn(ids_int, 2, simplify=F)
  pairs = c(pairs, lapply(seq_along(ids_int), function(i) c(i,i)))

  coefs = lapply(r, function(rr) {
    mem = as.environment(memTemplate)
    mem$r = rr

    unlist(lapply(pairs, function(p) {
      A = C = c(p[1], 0)
      B = D = c(p[2], 0)
      twoLocKin(A, B, C, D, mem, indent = ifelse(debug, 0, NA))
    }))
  })

  # Build result data frame
  labs = labels(x)
  idcols = do.call(rbind, pairs)
  idcols[] = labs[idcols]
  res = data.frame(id1 = idcols[,1], id2 = idcols[,2],
                   r = rep(r, each=length(pairs)),
                   phi2 = unlist(coefs),
                   stringsAsFactors = F)

  res
}


####################
# Internal functions
####################

twoLocKin = function(A, B, C, D, mem, indent = 0) {
  mem$i = mem$i + 1
  # Each A, B, C, D is a parent-child pair

  # If any founders -> 0
  if(A[1]*B[1]*C[1]*D[1] == 0) {print("Outside of pedigree! Can this happen?"); return(0)}

  # Sort: a >= b,c,d; c >= d; if(a==c) then b>=d
  plist = sortPairs(A,B,C,D)
  printMess(plist, indent)

  A = plist[[1]]; B = plist[[2]]; C = plist[[3]]; D = plist[[4]]
  a = A[1]; b = B[1]; c = C[1]; d = D[1]

  FF = mem$FIDX[a]
  MM = mem$MIDX[a]
  isFou = mem$isFounder[a]
  ANC = mem$anc
  k1 = mem$k1
  r = mem$r
  forceNonRec = a %in% mem$recomb$nr
  forceRec = a %in% mem$recomb$r

  ### The recursions (eqs. 7-11 in Thompson (1988))

  k2_recurse = function(A,B,C,D) {
    mem$irec = mem$irec + 1
    a = A[1]; b = B[1]; c = C[1]; d = D[1]

    if(a > b && a > c) { # eq. 7
      mem$eq7 = mem$eq7 + 1

      if(isFou)
        0
      else {
        t1 = twoLocKin(c(FF, a), B, C, D, mem, indent = indent + 2)
        t2 = twoLocKin(c(MM, a), B, C, D, mem, indent = indent + 2)
        0.5 * (t1 + t2)
      }
    }
    else if(a == b && a > c) { # eq. 8
      mem$eq8 = mem$eq8 + 1

      if(isFou)
        0.5 * k1[[c,d]]
      else {
        if(A[2] == B[2]) stop("this shouldn't happen")
        0.5 * (k1[[c,d]] + twoLocKin(c(FF, a), c(MM, a), C, D, mem, indent = indent + 2))
      }
    }
    else if(a > b && a == c) { # eqs. 9a & 9b
      if(A[2] == C[2]) { # eq. 9a
        mem$eq9a = mem$eq9a + 1

        if(isFou)
          0
        else if(forceNonRec) {
          t1 = twoLocKin(c(FF, a), B, c(FF, a), D, mem, indent = indent + 2)
          t2 = twoLocKin(c(MM, a), B, c(MM, a), D, mem, indent = indent + 2)

          0.5 * (1-r) * (t1 + t2)
        }
        else if(forceRec) {
          t3 = twoLocKin(c(FF, a), B, c(MM, a), D, mem, indent = indent + 2)
          t4 = twoLocKin(c(MM, a), B, c(FF, a), D, mem, indent = indent + 2)

          0.5 * r * (t3 + t4)
        }
        else {
          t1 = twoLocKin(c(FF, a), B, c(FF, a), D, mem, indent = indent + 2)
          t2 = twoLocKin(c(MM, a), B, c(MM, a), D, mem, indent = indent + 2)
          t3 = twoLocKin(c(FF, a), B, c(MM, a), D, mem, indent = indent + 2)
          t4 = twoLocKin(c(MM, a), B, c(FF, a), D, mem, indent = indent + 2)

          0.5 * ((1-r)*(t1 + t2) + r*(t3 + t4))
        }
      }
      else {  # eq. 9b
        mem$eq9b = mem$eq9b + 1
        if(isFou)
          0
        else {
          t1 = twoLocKin(c(FF, a), B, c(FF, a), D, mem, indent = indent + 2)
          t2 = twoLocKin(c(MM, a), B, c(MM, a), D, mem, indent = indent + 2)
          t3 = twoLocKin(c(FF, a), B, c(MM, a), D, mem, indent = indent + 2)
          t4 = twoLocKin(c(MM, a), B, c(FF, a), D, mem, indent = indent + 2)

          0.25 * (t1 + t2 + t3 + t4)
        }
      }
    }
    else if(a == b && a == c && a > d) { # eq. 10
      mem$eq10 = mem$eq10 + 1

      if(isFou)
        0
      else {
        t3 = twoLocKin(c(FF, a), c(MM, a), c(FF, a), D, mem, indent = indent + 2)
        t4 = twoLocKin(c(FF, a), c(MM, a), c(MM, a), D, mem, indent = indent + 2)
        s = 0.25 * (k1[[FF, d]] + k1[[MM, d]] + t3 + t4)

        if(forceNonRec)
          (1-r) * s
        else if(forceRec)
          r * s
        else
          s
      }
    }
    else if((a == b && a == c && a == d)) { # eq. 11
      mem$eq11 = mem$eq11 + 1
      R = .5*(r^2 + (1-r)^2)

      # This case needs further branching! A bit vague in the EAT-paper.
      type = sum(c(A[2], B[2]) %in% c(C[2], D[2]))

      if(type == 2) { # this is eq. 11 in the paper
        if(isFou) {
          if(forceNonRec && forceRec) 0
          else if(forceNonRec) .5 * (1-r)^2
          else if(forceRec) .5 * r^2
          else R
        }
        else {
          if(forceNonRec && forceRec) r*(1-r)*k1[[MM, FF]] # without the factor two (either R-NR or NR-R)
          else {
            t4 = twoLocKin(c(MM, a), c(FF, a), c(MM, a), c(FF, a), mem, indent = indent + 2)
            if(forceNonRec) .5 * (1-r)^2 * (1 + t4)
            else if(forceRec) .5 * r^2 * (1 + t4)
            else 2*r*(1-r)*k1[[MM, FF]] + R*(1 + t4)
          }
        }
      }
      else if(type == 1) { # a.k.a. k2(J(A1, A2), L(A1, A3))
        if(forceNonRec | forceRec)
          stop2("Special case: Not implemented.")
        if(isFou)
          1/4
        else {
          t1 = twoLocKin(c(FF, a), B, c(FF, a), D, mem, indent = indent + 2)
          t2 = twoLocKin(c(MM, a), B, c(MM, a), D, mem, indent = indent + 2)
          t3 = twoLocKin(c(FF, a), B, c(MM, a), D, mem, indent = indent + 2)
          t4 = twoLocKin(c(MM, a), B, c(FF, a), D, mem, indent = indent + 2)

          0.5 * ((1-r)*(t1 + t2) + r*(t3 + t4))
        }
      }
      else if(type == 0) { # a.k.a. k2(J(A1, A2), L(A3, A4))
        stop2("Special case k2(J(A1, A2), L(A3, A4)): Not implemented yet")
      }
    }
  }

  # If no common ancestors, return 0
  if(!(ANC[a,b] && ANC[c,d])) {
    if(!is.na(indent)) message(strrep(" ", indent), 0)
    return(0)
  }

  # Lookup in array; compute if necessary.
  mem$ilook = mem$ilook + 1
  res = mem$k2[[toString(plist)]]
  if(is.null(res))
    res = mem$k2[[toString(plist)]] = k2_recurse(A, B, C, D)

  if(!is.na(indent)) message(strrep(" ", indent), res)
  res
}

sortPairs = function(A,B,C,D) {
  # Sort: a >= b,c,d; c >= d; if(a==c) then b>=d
  plist = list(A,B,C,D)

  if(plist[[1]][1] < plist[[2]][1])
    plist[1:2] = plist[2:1]
  if(plist[[3]][1] < plist[[4]][1])
    plist[3:4] = plist[4:3]
  if(plist[[1]][1] < plist[[3]][1] ||
     (plist[[1]][1] == plist[[3]][1] && plist[[2]][1] < plist[[4]][1]))
    plist[] = plist[c(3,4,1,2)]
  plist
}

printMess = function(plist, indent) {
  if(is.na(indent)) return()
  pp = sapply(plist, function(p) paste(p, collapse="-"))
  message(sprintf("%sA = %s, B = %s, C = %s, D = %s",
          strrep(" ", indent), pp[1], pp[2], pp[3], pp[4]))
}

initialiseTwoLocusMemo = function(ped, r, recomb = NULL, chromType = "autosomal") {
  # Create memory storage
  mem = new.env()

  # Start timing
  st = Sys.time()

  mem$FIDX = ped$FIDX
  mem$MIDX = ped$MIDX
  mem$SEX = ped$SEX
  mem$r = r

  # Conditions on recombinant/non-recombinant gametes
  if(is.null(recomb))
    recomb = list(r = NULL, nr = NULL)
  mem$recomb = recomb

  # Logical matrix showing who has a common ancestor within the pedigree.
  mem$anc = has_common_ancestor(ped)

  # Compute kinship matrix directly
  mem$k1 = switch(chromType, autosomal = kinship(ped), x = kinshipX(ped))

  # Storage for two-locus kinship values
  mem$k2 = list()

  # For quick look-up:
  FOU = founders(ped, internal = T)
  isFounder = rep(FALSE, pedsize(ped))
  isFounder[FOU] = TRUE
  mem$isFounder = isFounder

  # Counters
  mem$i = mem$ilook = mem$irec = 0
  mem$eq7 = mem$eq8 = mem$eq9a = mem$eq9b = mem$eq10 = mem$eq11 = 0

  # Start time
  mem$st = st
  mem$initTime = Sys.time() - st

  mem
}

