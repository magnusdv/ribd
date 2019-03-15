#' Two-locus kinship coefficients
#'
#' Computes the two-locus kinship coefficient for two pedigree members, at two
#' loci separated by a given recombination rate. WORK IN PROGRESS
#'
#' The implementation is based on the recursive algorithm described by Thompson
#' (1988).
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param r A number in the interval \eqn{[0, 0.5]}; the recombination rate between the two loci.
#' @param verbose A logical
#'
#'
#' @references E. A. Thompson (1988). _Two-locus and Three-locus Gene Identity
#'   by Descent in Pedigrees_ IMA Journal of Mathematics Applied in Medicine &
#'   Biology, vol. 5.
#'
#' @examples
#' # Full sibs
#' x = nuclearPed(2)
#'
#' k_0 = twoLocusKinship(x, ids = 3:4, r = 0)
#' k_0.5 = twoLocusKinship(x, ids = 3:4, r = 0.5)
#'
#' stopifnot(k_0 == 1/4, k_0.5 == 1/16)
#'
#' #######################################
#' # Reproducing Fig. 3 in Thompson (1988)
#' # Note that in the article, curve (a) is wrong.
#' # (See published Erratum.)
#' #######################################
#'
#' # Pedigrees (a) - (d)
#' peds = list(
#'   a = list(ped = linearPed(3), ids = c(1,7)),
#'   b = list(ped = halfCousinPed(0, 1), ids = c(3,7)),
#'   c = list(ped = cousinPed(1), ids = c(5,8)),
#'   d = list(ped = doubleCousins(1, 1, half1 = TRUE, half2 = TRUE),
#'            ids = c(5,9)))
#'
#' # Recombination values
#' rseq = seq(0, 0.5, length = 20)
#'
#' # Compute two-locus kinship coefficients
#' kvals = sapply(peds, function(x)
#'   sapply(rseq, function(r) twoLocusKinship(x$ped, x$ids, r)))
#'
#' # Plot
#' matplot(rseq, kvals, type="l")
#' legend("topright", letters[1:4], col = 1:4, lty = 1:4)
#'
#' @export
twoLocusKinship = function(x, ids, r, verbose = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  mem = initialiseTwoLocusMemo(x, r)

  ids_int = internalID(x, ids)
  id1 = ids_int[1]
  id2 = ids_int[2]
  A = C = c(id1, 0)
  B = D = c(id2, 0)
  res = twoLocKin(A, B, C, D, mem, indent = ifelse(verbose, 0, NA))

  # Print info
  if(verbose) {
    initsecs = sprintf("%.2f", mem$initTime)
    totsecs = sprintf("%.1f", Sys.time()-mem$st)
    print(glue::glue("
                     Calls = {mem$i}
                     Lookups = {mem$ilook}
                     Recursions: {mem$irec}
                     Total time used: {totsecs} seconds"))
  }

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

  ### The recursions (eqs. 7-11 in Thompson (1988))

  k2_recurse = function(A,B,C,D) {
    mem$irec = mem$irec + 1
    a = A[1]; b = B[1]; c = C[1]; d = D[1]

    if(a > b && a > c) { # eq. 7
      if(isFou)
        0
      else {
        t1 = twoLocKin(c(FF, a), B, C, D, mem, indent = indent + 2)
        t2 = twoLocKin(c(MM, a), B, C, D, mem, indent = indent + 2)
        0.5 * (t1 + t2)
      }
    }
    else if(a == b && a > c) { # eq. 8
      if(isFou)
        0.5 * k1[[c,d]]
      else {
        if(A[2] == B[2]) stop("this shouldn't happen")
        0.5 * (k1[[c,d]] + twoLocKin(c(FF, a), c(MM, a), C, D, mem, indent = indent + 2))
      }
    }
    else if(a > b && a == c) { # eqs. 9a, 9b
      if(isFou)
        0
      else {
        t1 = twoLocKin(c(FF, a), B, c(FF, a), D, mem, indent = indent + 2)
        t2 = twoLocKin(c(MM, a), B, c(MM, a), D, mem, indent = indent + 2)
        t3 = twoLocKin(c(FF, a), B, c(MM, a), D, mem, indent = indent + 2)
        t4 = twoLocKin(c(MM, a), B, c(FF, a), D, mem, indent = indent + 2)

        if(A[2] == C[2]) { # eq. 9a
          0.5 * ((1-r)*(t1 + t2) + r*(t3 + t4))
        }
        else { # eq. 9a
          0.25 * (t1 + t2 + t3 + t4)
        }
      }
    }
    else if(a == b && a == c && a > d) { # eq. 10
      if(isFou)
        0
      else {
        t3 = twoLocKin(c(FF, a), c(MM, a), c(FF, a), D, mem, indent = indent + 2)
        t4 = twoLocKin(c(FF, a), c(MM, a), c(MM, a), D, mem, indent = indent + 2)
        0.25 * (k1[[FF, d]] + k1[[MM, d]] + t3 + t4)
      }
    }
    else if((a == b && a == c && a == d)) { # eq. 11
      R = .5*(r^2 + (1-r)^2)
      if(isFou)
        R
      else {
        t4 = twoLocKin(c(MM, a), c(FF, a), c(MM, a), c(FF, a), mem, indent = indent + 2)
        2*r*(1-r)*k1[[MM, FF]] + R*(1 + t4)
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

initialiseTwoLocusMemo = function(ped, r, chromType = "autosomal") {
  # Create memory storage
  mem = new.env()

  # Start timing
  st = Sys.time()

  mem$FIDX = ped$FIDX
  mem$MIDX = ped$MIDX
  mem$SEX = ped$SEX
  mem$r = r

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

  # Start time
  mem$st = st
  mem$initTime = Sys.time() - st

  mem
}

test = function(x) {
  if(!is.ped(x)) return()
  k1 = function(x, ids) kinship2_kinship(x, ids)
  k2_0 = function(x,ids) twoLocusKinship(x,ids,r=0, verbose = F)
  k2_05 = function(x,ids) twoLocusKinship(x,ids,r=0.5, verbose = F)
  for(ids in combn(labels(x), 2, simplify = F)) {
    if(k2_0(x,ids) != k1(x,ids)) {
      cat("Error r=0!", ids, "\n")
      #x11(); plot(x, col=list(red=ids))
      #return(list(x, ids))
    }
    if(k2_05(x,ids) != k1(x,ids)^2) {
      cat("Error r=0.5!", ids, "\n")
      #x11(); plot(x, col=list(red=ids))
      return(list(x, ids))
    }
  }
}

