#' Generalised kinship coefficients
#'
#' Compute generalised single-locus kinship coefficients, using the algorithm by
#' Weeks & Lange (1992).
#' @param x A `ped` object.
#' @param pattern A `kinPattern` object.
#' @param mem An `environment` (for internal use).
#' @param verbose A logical.
#' @param debug A logical.
#'
#' @return A single probability.
#'
#' @references Weeks and Lange (1992):
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
#' y = quadHalfFirstCousins()
#' ids = 9:10
#' stopifnot(IBD_from_gk(y, ids) == kappaIBD(y, ids))
#'
#' y = nuclearPed(2); ids = 3:4
#' stopifnot(IBD_from_gk(y, ids) == kappaIBD(y, ids))
#'
#' #### Triple/quad kinship (compare with karigl)
#' x = fullSibMating(1)
#' ids = c(1,5,6)
#' stopifnot(generalisedKinship(x, list(ids)) == generalisedKinship3(x, ids))
#' ids = c(1,5,6,5)
#' stopifnot(generalisedKinship(x, list(ids)) == generalisedKinship4(x, ids))
#'
#' @export
generalisedKinship = function(x, pattern, mem = NULL, verbose = F, debug = F) {

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  x = foundersFirst(x)

  if(!inherits(pattern, "kinPattern") && is.list(pattern))
    pattern = kinPattern(x, pattern)

  # Setup memoisation
  if(is.null(mem))
    mem = initialiseGKMemo(x, counters = c("i", "itriv", "iimp", "ifound", "ilook", "irec"))

  res = genKin(pattern, mem, indent = ifelse(debug, 0, NA))

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
  uniqVec = unlist(uniqList, use.names = F)

    # Boundary 1: Anyone in >2 groups?
  if(any(tabulate(uniqVec) > 2)) {
    mem$iimp = mem$iimp + 1
    return(printAndReturn(0, indent, comment = " (B1)"))
  }


  # Boundary 2: Any group with 2 unrelated indivs?
  k1 = mem$k1
  for(s in uniqList[lengths(uniqList) > 1]) {
    pairs_mat = comb2(s, vec = T)
    if(any(k1[pairs_mat] == 0)) {
      mem$iimp = mem$iimp + 1
      return(printAndReturn(0, indent, comment = " (B2)"))
    }
  }


  if(all(mem$isFounder[unlist(uniqList)])) {
    mem$ifound = mem$ifound + 1
    m1 = sum(lengths(kp))
    m2 = length(unique.default(uniqVec))
    val = 1/2^(m1 - m2)
    return(printAndReturn(val, indent, comment = " (B3)"))
  }

  kp = sort_kinPattern(kp)

  # Lookup in array; compute if necessary.
  kinStr = toString.default(kp)
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

  s = sum(kp[[1]] == pivot)
  t = if(length(kp) > 1) sum(kp[[2]] == pivot) else 0

  # Recurse!
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
  if(t > 0) {
    A2 = .5^(s + t)
    res =
      A2 * recu(kinReplace1(kp, id = pivot, rep1 = fa, rep2 = mo)) +
      A2 * recu(kinReplace1(kp, id = pivot, rep1 = mo, rep2 = fa))
  }

  mem$PHI[[kinStr]] = res
  return(printAndReturn(res, indent))
}



# Constructor for S3 class "kinPattern"
newKinPattern = function(pattern, labels) {
  if(!is.list(pattern) || !is.numeric(unlist(pattern)) || !is.character(labels))
    stop2("Invalid input")

  pattern = lapply(pattern, as.integer)
  structure(pattern, labels = labels, class = "kinPattern")
}


# Helper function for creating kinPattern objects
#' @export
kinPattern = function(x, pattern, internal = F) {
  if(!is.ped(x))
    stop2("First argument must be a `ped` object")

  if(!internal) {
    pattern = lapply(pattern, internalID, x = x)
  }
  else {
    v = unlist(pattern)
    if(!is.numeric(v))
      stop2("Non-numeric entries found in `pattern`; this is illegal when `internal = T`")
    if(length(err <- setdiff(v, 1:pedsize(x))) > 0)
      stop2("Illegal entry in kinPattern() when `internal=T`: ", err)
  }

  newKinPattern(pattern, labels(x))
}


print.kinPattern = function(x, ...,  indent = 0) {
  if(is.na(indent))
    return()
  labs = attr(x, "labels")

  grps = vapply(x, function(g) sprintf("(%s)", paste0(labs[g], collapse = ",")), FUN.VALUE = "")
  str = paste0(grps, collapse = ",")

  cat(strrep(" ", indent), str, "\n", sep = "")
}


sort_kinPattern = function(x) {

  # Sort each group
  x[] = lapply(x, function(g) sort.int(g, decreasing = T, method = "shell"))

  # Max number of elements
  L = max(lengths(x))

  # Max ID
  M = max(unlist(x, use.names = F))

  # Convert to number, used for sorting
  b = vapply(x, function(g) sum(g * ((M+1)^(seq(L, by = -1, length = length(g))))), FUN.VALUE = 1)

  x[] = x[order(b, decreasing = T, method = "shell")]
  x
}


kinReduce1 = function(x) {
  # Remove empty groups
  x[lengths(x) == 0] = NULL
  x
}


kinReplace1 = function(x, id, rep1, rep2 = NULL) {
  g1 = x[[1]]
  x[[1]] = c(rep1, g1[g1 != id])

  if(!is.null(rep2)) {
    g2 = x[[2]]
    x[[2]] = c(rep2, g2[g2 != id])
  }
  x
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
  # TODO: is this neccessary? (since we also include k1 below)
  mem$anc = hasCommonAncestor(ped)

  # Compute kinship matrix directly
  mem$k1 = switch(chromType, autosomal = kinship(ped), x = kinshipX(ped))

  # Storage for result values
  mem$PHI = list()

  # For quick look-up:
  FOU = founders(ped, internal = T)
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

