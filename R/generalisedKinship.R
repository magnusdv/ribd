#' Generalised kinship coefficients
#'
#' Compute generalised single-locus kinship coefficients, as defined by Weeks &
#' Lange (1988), and also the further generalisation by Lange & Sinsheimer
#' (1992) allowing to specify if the alleles involved are paternally or
#' maternally inherited.
#'
#' @param x A `ped` object.
#' @param pattern A `kinPattern` object.
#' @param mem An `environment` (for internal use).
#' @param verbose A logical.
#' @param debug A logical.
#'
#' @return A single probability.
#'
#' @references Weeks and Lange. _The Affected-Pedigree-Member Method of Linkage
#'   Analysis_. Am. J. Hum. Genet. 42:315-326, 1988.
#'
#' @examples
#' x = nuclearPed(3)
#' kp = kinPattern(x, list(c(1,1,1)))
#' generalisedKinship(x, kp)
#'
#'
#' ##### IBD coefficients via generalised kinship ###
#' #(Clearly not the simplest way; serves as a check)
#' IBD_from_gk = function(x, ids) {
#'   fa1 = father(x, ids[1])
#'   fa2 = father(x, ids[2])
#'   mo1 = mother(x, ids[1])
#'   mo2 = mother(x, ids[2])
#'   GK = function(...) generalisedKinship(x, list(...))
#'
#'   k0 = GK(fa1, fa2, mo1, mo2)
#'   k1 = GK(c(fa1, fa2), mo1, mo2) + GK(c(fa1, mo2), fa2, mo1) +
#'        GK(c(mo1, fa2), fa1, mo2) + GK(c(mo1, mo2), fa1, fa2)
#'   k2 = GK(c(fa1, fa2), c(mo1, mo2)) + GK(c(fa1, mo2), c(mo1, fa2))
#'   c(k0, k1, k2)
#' }
#'
#' y1 = nuclearPed(2); ids = 3:4
#' stopifnot(IBD_from_gk(y1, ids) == kappaIBD(y1, ids))
#'
#' y2 = quadHalfFirstCousins()
#' ids = 9:10
#' stopifnot(IBD_from_gk(y2, ids) == kappaIBD(y2, ids))
#'
#' #### Triple/quad kinship (compare with karigl)
#' x = fullSibMating(1)
#' ids = c(1,5,6)
#' stopifnot(generalisedKinship(x, list(ids)) == generalisedKinship3(x, ids))
#' ids = c(1,5,6,5)
#' stopifnot(generalisedKinship(x, list(ids)) == generalisedKinship4(x, ids))
#'
#' @export
generalisedKinship = function(x, pattern, mem = NULL, verbose = FALSE, debug = FALSE) {

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  x = foundersFirst(x)

  # Need both original and internal version
  if(inherits(pattern, "kinPattern")) {
    kp = pattern
    pattern = as.list.kinPattern(kp)
  }
  else {
    pattern = unname(pattern) # remove list names if present
    kp = kinPattern(x, pattern)
  }

  det = isDetailed(kp)

  # If founders, adjust pedigree by adding parents
  if(det) {
    for(id in intersect(unlist(pattern), founders(x)))
      x = addParents(x, id, verbose = FALSE)

    x = foundersFirst(x)
    kp = kinPattern(x, pattern)
  }

  # Setup memoisation
  if(is.null(mem))
    mem = initialiseGKMemo(x, counters = c("i", "itriv", "iimp", "ifound", "ilook", "irec"))

  gkFUN = if(isDetailed(kp)) genKinDetailed else genKin
  res = gkFUN(kp, mem, indent = ifelse(debug, 0, NA))

  if(verbose)
    printCounts2(mem)

  res
}


genKin = function(kp, mem, indent = 0) {
  mem$i = mem$i + 1

  print.kinPattern(kp, indent = indent)
  kp = kinReduce1(kp)

  if(length(kp) == 0 || (length(kp) == 1 && length(kp[[1]]) == 1)) {
    mem$itriv = mem$itriv + 1
    return(printAndReturn(1, indent, comment = " (empty)"))
  }

  uniqList = lapply(kp, unique.default)
  uniqVec = unlist(uniqList, use.names = FALSE)

    # Boundary 1: Anyone in >2 groups?
  if(any(tabulate(uniqVec) > 2)) {
    mem$iimp = mem$iimp + 1
    return(printAndReturn(0, indent, comment = " (B1)"))
  }

  # Boundary 2: Any group with 2 unrelated indivs?
  k1 = mem$k1
  for(s in uniqList[lengths(uniqList) > 1]) {
    pairs_mat = comb2(s, vec = TRUE)
    if(any(k1[pairs_mat] == 0)) {
      mem$iimp = mem$iimp + 1
      return(printAndReturn(0, indent, comment = " (B2)"))
    }
  }

  # Boundary 3: Only founders in each group?
  if(all(mem$isFounder[unlist(uniqList)])) {
    mem$ifound = mem$ifound + 1
    m1 = sum(lengths(kp))
    m2 = length(unique.default(uniqVec))
    val = 1/2^(m1 - m2)
    return(printAndReturn(val, indent, comment = " (B3)"))
  }

  kp = sort_kinPattern(kp)

  # Lookup in array; compute if necessary.
  kinStr = format.kinPattern(kp)
  val = mem$PHI[[kinStr]]

  if(!is.null(val)) {
    mem$ilook = mem$ilook + 1
    mem$PHI[[kinStr]] = val
    return(printAndReturn(val, indent, comment = " (lookup)"))
  }

  mem$irec = mem$irec + 1

  # Wrapper of genKin to save typing in the recursions
  recu = function(k) genKin(k, mem, indent = indent +2)

  pivot = kp[[1]][1]
  fa = mem$FIDX[pivot]
  mo = mem$MIDX[pivot]

  # Number of times the pivot occurs in first block (s) and second block (t)
  s = sum(kp[[1]] == pivot)
  t = if(length(kp) > 1) sum(kp[[2]] == pivot) else 0

  # Recurrence rules 1 (s = 1) and 2
  if(t == 0) {
    A1 = .5^s
    res =
      A1 * recu(kinReplace1(kp, id = pivot, rep1 = fa)) +
      A1 * recu(kinReplace1(kp, id = pivot, rep1 = mo))

    if(s > 1) {
      B1 = (1 - 2 * (.5)^s)
      res = res +
        B1 * recu(kinReplace1(kp, id = pivot, rep1 = c(fa, mo)))
    }
  }

  # Recurrence rule 3
  if(t > 0) {
    A2 = .5^(s + t)
    res =
      A2 * recu(kinReplace1(kp, id = pivot, rep1 = fa, rep2 = mo)) +
      A2 * recu(kinReplace1(kp, id = pivot, rep1 = mo, rep2 = fa))
  }

  mem$PHI[[kinStr]] = res
  return(printAndReturn(res, indent))
}



#' Generalised kinship pattern
#'
#' @param x A `ped` object
#' @param pattern A list of vectors of ID labels.
#' @param internal A logical
#'
#' @return An object of class `kinPattern`.
#'
#' @examples
#' kinPattern(nuclearPed(2), list(1, 3:4))
#'
#' @export
kinPattern = function(x, pattern, internal = FALSE) {
  if(!is.ped(x))
    stop2("First argument must be a `ped` object")

  v = unlist(pattern)
  detailed = !is.null(allnms <- names(v))
  if(detailed && !all(allnms %in% c("", "p", "m")))
    stop2("Names (to indicate deterministic sampling) must be 'p' or 'm': ",
          setdiff(allnms, c("", "p", "m")))

  if(!internal) {
    pattern = lapply(pattern, function(v) {
      w = internalID(x, v)
      if(detailed && !is.null(nms <- names(v)) && !all(nms == ""))
        names(w) = names(v)
      w
    })
  }
  else {
    if(!is.numeric(v))
      stop2("Non-numeric entries found in `pattern`; this is illegal when `internal = TRUE`")
    if(length(err <- setdiff(v, 1:pedsize(x))) > 0)
      stop2("Illegal entry in kinPattern() when `internal = TRUE`: ", err)
    pattern = lapply(pattern, `mode<-`, "integer")
  }

  newKinPattern(pattern, labels(x))
}

# Constructor for S3 class "kinPattern"
newKinPattern = function(pattern, labels) {
  if(!is.list(pattern) || !is.numeric(unlist(pattern)) || !is.character(labels))
    stop2("Invalid input")

  pattern = lapply(pattern, `mode<-`, "integer")

  nms = names(unlist(pattern))
  detailed = !is.null(nms)

  # If pattern has all-empty names, convert to NULL
  if(detailed && all(nms == "")) {
    pattern = lapply(pattern, `names<-`, NULL)
    detailed = FALSE
  }

  structure(pattern, labels = labels, detailed = detailed, class = "kinPattern")
}

format.kinPattern = function(x, ...) {
  labs = attr(x, "labels")
  detailed = isDetailed(x)

  grps = vapply(x, function(g) {
    lb = labs[g]
    if(detailed) {
      idx = names(g) != ""
      lb[idx] = paste(lb[idx], names(g)[idx], sep = ":")
    }
    sprintf("(%s)", paste0(lb, collapse = ","))
  }, FUN.VALUE = "")

  paste0(grps, collapse = ",")
}

print.kinPattern = function(x, ..., indent = 0) {
  if(is.na(indent))
    return()

  str = format.kinPattern(x)
  cat(strrep(" ", indent), str, "\n", sep = "")
}

as.list.kinPattern = function(x) {
  labs = attr(x, "labels")
  det = isDetailed(x)

  lapply(x, function(v) {
      w = labs[v]
      if(det) names(w) = names(v)
      w
    })
}

sort_kinPattern = function(x) {

  det = isDetailed(x)

  # Function for sorting a single group
  sortFUN = function(g) {
      if(!det || is.null(nms <- names(g)))
        sort.int(g, decreasing = TRUE, method = "shell")
      else
        g[order(g, nms, decreasing = TRUE, method = "shell")]
  }

  # Sort each group
  x[] = lapply(x, sortFUN)

  # Max number of elements
  L = max(lengths(x))

  # Max ID
  M = max(unlist(x, use.names = FALSE))

  # Convert to number, used for sorting
  b = vapply(x, function(g) sum(g * ((M+1)^(seq(L, by = -1, length = length(g))))), FUN.VALUE = 1)

  # Overall group order
  x[] = x[order(b, decreasing = TRUE, method = "shell")]

  x
}

isDetailed = function(kp) {
  det = attr(kp, "detailed")
  !is.null(det) && det
}

kinReduce1 = function(kp) {
  # Remove empty groups
  kp[lengths(kp) == 0] = NULL
  kp
}

kinReduce1Det = function(kp) {

  stillDetailed = FALSE

  # Remove deterministic repeats in each block
  for(i in seq_along(kp)) {
    g = kp[[i]]
    nms = names(g)
    if(is.null(nms))
      next
    det = nms != ""
    if(!any(det)) {
      names(kp[[i]]) = NULL
      next
    }

    stillDetailed = TRUE

    dups = duplicated.default(paste(g[det], nms[det]))
    if(any(dups)) {
      remov = which(det)[dups]
      kp[[i]] = g[-remov]
    }
  }

  if(!stillDetailed)
    attr(kp, "detailed") = FALSE

  kinReduce1(kp)
}


kinReplace1 = function(kp, id, rep1, rep2 = NULL) {
  g1 = kp[[1]]
  kp[[1]] = c(rep1, g1[g1 != id])

  if(!is.null(rep2)) {
    g2 = kp[[2]]
    kp[[2]] = c(rep2, g2[g2 != id])
  }
  kp
}


# Initialise memoisation used in twoLocusKinship and twoLocusGeneralisedKinship
initialiseGKMemo = function(ped, chromType = "autosomal", counters = NULL) {
  if(chromType != "autosomal")
    stop2("Only `chromType = autosomal` is implemented at the moment")

  # Create memory storage
  mem = new.env()

  # Start timing
  st = Sys.time()

  mem$FIDX = ped$FIDX
  mem$MIDX = ped$MIDX
  mem$SEX = ped$SEX

  # Logical matrix showing who has a common ancestor within the pedigree.
  # TODO: is this necessary? (since we also include k1 below)
  # mem$anc = hasCommonAncestor(ped)

  # Compute kinship matrix directly
  mem$k1 = kinship(ped, Xchrom = chromType == "x")

  # Storage for result values
  mem$PHI = list()

  # For quick look-up:
  FOU = founders(ped, internal = TRUE)
  isFounder = rep(FALSE, pedsize(ped))
  isFounder[FOU] = TRUE
  mem$isFounder = isFounder

  ### Founder inbreeding
  #TODO

  # Counters
  for(cou in counters)
    assign(cou, 0, envir = mem)

  # Start time
  mem$st = st
  mem$initTime = Sys.time() - st

  mem
}



genKinDetailed = function(kp, mem, indent = 0) {
  mem$i = mem$i + 1

  print.kinPattern(kp, indent = indent)
  kp = kinReduce1Det(kp)

  if(length(kp) == 0 || (length(kp) == 1 && length(kp[[1]]) == 1)) {
    mem$itriv = mem$itriv + 1
    return(printAndReturn(1, indent, comment = " (empty)"))
  }

  # If no longer detailed, use condensed function
  if(!isDetailed(kp))
    return(genKin(kp, mem, indent = indent))

  uniqList = lapply(kp, unique.default)
  uniqVec = unlist(uniqList, use.names = FALSE)

  ### Boundary functions: same as in the condensed case

  # Boundary 1: Anyone in >2 groups?
  counts = tabulate(uniqVec, nbins = max(uniqVec))
  if(any(counts > 2)) {
    mem$iimp = mem$iimp + 1
    return(printAndReturn(0, indent, comment = " (B1)"))
  }

  # Boundary 2: Any group with 2 unrelated indivs?
  k1 = mem$k1
  for(s in uniqList[lengths(uniqList) > 1]) {
    pairs_mat = comb2(s, vec = TRUE)
    if(any(k1[pairs_mat] == 0)) {
      mem$iimp = mem$iimp + 1
      return(printAndReturn(0, indent, comment = " (B2)"))
    }
  }

  # Boundary 3: Only founders in each group?
  if(all(mem$isFounder[uniqVec])) {
    mem$ifound = mem$ifound + 1
    m1 = sum(lengths(kp))
    m2 = length(unique.default(uniqVec))
    val = 1/2^(m1 - m2)
    return(printAndReturn(val, indent, comment = " (B3)"))
  }

  # Restriction D1: Anyone with paternal (resp. maternal) allele in multiple blocks?
  # (Recall: kp is reduced, so nobody has repeated det-sampling *within* a block)
  kpFlat = unlist(kp, use.names = TRUE)
  nmsFlat = names(kpFlat)

  if(anyDuplicated.default(kpFlat[nmsFlat == "p"]) || anyDuplicated.default(kpFlat[nmsFlat == "m"])) {
    mem$iimp = mem$iimp + 1
    return(printAndReturn(0, indent, comment = " (D1)"))
  }

  # Restriction D2: Anyone with mat&pat in same block AND random in other block?
  kpDet = kpFlat[nmsFlat != ""]
  both = kpDet[duplicated.default(kpDet)] # all with both p and m

  grp = rep(seq_along(kp), lengths(kp)) # group index vector

  for(i in both[counts[both] == 2]) { # all with both AND in two groups
    grpDet = grp[kpFlat == i & nmsFlat != ""] # groups where i is deterministic
    if(grpDet[1] == grpDet[2]) {  # check if m&p in same group
      mem$iimp = mem$iimp + 1
      return(printAndReturn(0, indent, comment = " (D2)"))
    }
  }

  kp = sort_kinPattern(kp)

  # Lookup in array; compute if necessary.
  kinStr = format.kinPattern(kp)
  val = mem$PHI[[kinStr]]

  if(!is.null(val)) {
    mem$ilook = mem$ilook + 1
    mem$PHI[[kinStr]] = val
    return(printAndReturn(val, indent, comment = " (lookup)"))
  }

  mem$irec = mem$irec + 1

  # Wrapper of genKin to save typing in the recursions
  recu = function(k) genKinDetailed(k, mem, indent = indent +2)

  g1 = kp[[1]]
  g2 = if(length(kp) > 1) kp[[2]] else NULL
  nms1 = names(g1)
  nms2 = names(g2)

  # Pivot individual (after sorting: highest)
  pivot = g1[1]
  fa = mem$FIDX[pivot]
  mo = mem$MIDX[pivot]

  # Number of times the pivot occurs in first block (s) and second block (t)
  s = sum(g1 == pivot)
  t = sum(g2 == pivot)

  # Is deterministic sampling involved in first group?
  det1 = nms1[g1 == pivot & nms1 != ""]  # NULL or vector of length 0, 1 or 2
  s = s - length(det1)

  # Recurrence rule 1a,b,c
  if(t == 0) {

    A1 = .5^s

    if(length(det1) == 0) { # Rule 1a
      res =
        A1 * recu(kinReplace1(kp, id = pivot, rep1 = fa)) +
        A1 * recu(kinReplace1(kp, id = pivot, rep1 = mo)) +
        if(s > 1) (1 - 2*A1) * recu(kinReplace1(kp, id = pivot, rep1 = c(fa, mo))) else 0
    }
    else if(length(det1) == 1) { # Rule 1b
      res =
        A1 * recu(kinReplace1(kp, id = pivot, rep1 = switch(det1, p = fa, m = mo))) +
        if(s > 0) (1 - A1) * recu(kinReplace1(kp, id = pivot, rep1 = c(fa, mo))) else 0
    }
    else if(length(det1) == 2) { # Rule 1c
      res = recu(kinReplace1(kp, id = pivot, rep1 = c(fa, mo)))
    }
  }
  else {  # t > 0

    # Deterministic sampling in second group?
    det2 = nms2[g2 == pivot & nms2 != ""]  # NULL or vector of length 0, 1 or 2
    t = t - length(det2)

    A2 = .5^(s + t)

    if(length(det1) == 0 && length(det2) == 0) { # Rule 2a
    res =
      A2 * recu(kinReplace1(kp, id = pivot, rep1 = fa, rep2 = mo)) +
      A2 * recu(kinReplace1(kp, id = pivot, rep1 = mo, rep2 = fa))
    }
    else {  # Rules 2b and 2c
      rep1 = if(length(det1) == 1) switch(det1, p = fa, m = mo) else switch(det2, p = mo, m = fa)
      rep2 = setdiff(c(fa, mo), rep1)
      res = A2 * recu(kinReplace1(kp, id = pivot, rep1 = rep1, rep2 = rep2))
    }
  }

  mem$PHI[[kinStr]] = res
  return(printAndReturn(res, indent))
}

