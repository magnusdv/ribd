# Constructor for S3 class "kin2L"
newKin2L= function(locus1, locus2, labels) {
  if(!is.list(locus1) || !is.list(locus2) || !is.character(labels))
    stop2("Invalid input to `newKinshipClass()`")

  for(g in c(locus1, locus2)) {
    if(!setequal(names(g), c("from", "to")))
      stop2("Group names must be 'from' and 'to'")
    if(!is.integer(g$from) || !is.integer(g$to))
      stop2("Group entries must be integer vectors, but are: ", class(g$from), ", ", class(g$to))
  }

  structure(list(locus1 = locus1, locus2 = locus2), labels = labels, class = "kin2L")
}

# Helper function for creating kin2L objects
kin2L = function(x, locus1, locus2, internal = F) {
  if(!is.ped(x))
    stop2("First argument must be a `ped` object")

  # If character input, convert to list. Example: "5>9 = 7>10, 6>9, 8>10"
  if(is.character(locus1)) locus1 = char2kinList(locus1)
  if(is.character(locus2)) locus2 = char2kinList(locus2)

  if(!internal) {
    locus1 = kinList2internal(x, locus1)
    locus2 = kinList2internal(x, locus2)
  }
  else { # ensure integer entries
    locus1 = rapply(locus1, as.integer, how = "replace")
    locus2 = rapply(locus2, as.integer, how = "replace")
  }

  obj = newKin2L(locus1, locus2, labels(x))
  validateKin2L(obj, ped = x)
}


# Validate kin2L objects
validateKin2L = function(x, ped = NULL) {

  printAndStop = function(g, err) {print(g); stop2(err)}

  for(i in 1:2) {
    loc = x[[i]]
    for(g in loc) {
      if(!is.list(g) || length(g) != 2)
        printAndStop(g, "Allele group is not a list of length 2")
      if(!setequal(names(g), c("from", "to")))
        printAndStop(g, "Allele group must have names 'from' and 'to'")
      if(!is.integer(g$from) || !is.integer(g$to) || length(g$from) != length(g$to))
        printAndStop(g, "Entries 'from' and 'to' must be numeric vectors (internal labels) of the same length")

      # # Check that meioses are compatible with pedigree
      # target.in.ped = g$to > 0
      # if(!is.null(ped) && length(target.in.ped) > 0) {
      #   labs = labels(ped)
      #   ch = g$to[target.in.ped]
      #   par = g$from[target.in.ped]
      #   par.ok = par == ped$FIDX[ch] | par == ped$MIDX[ch]
      #   if(!all(par.ok))
      #     printAndStop(x, toString(sprintf("%s is not a child of %s", labs[ch], labs[par])))
      # }
    }
  }
  x
}

print.kin2L = function(x, ...,  indent = 0) {
  if(is.na(indent)) return()
  labs = attr(x, "labels")

  mess1 = kinList2char(x$locus1, labs = labs)
  mess2 = kinList2char(x$locus2, labs = labs)

  message(sprintf("%s[%s ::: %s]", strrep(" ", indent), mess1, mess2))
}


sort_kin2L = function(x) {
  # Auxiliary function for sorting a single group
  sortGroup = function(g) {
    ord = order(g$from, -g$to, decreasing = T)
    list(from = g$from[ord], to = g$to[ord])
  }

  # Find highest ID
  x1 = x[[1]]
  x2 = x[[2]]
  pivot = max(unlist(lapply(c(x1, x2), function(g) g$from), use.names = F))

  # Which groups are pivot groups?
  piv1 = vapply(x1, function(g) pivot %in% g$from, FUN.VALUE = F)
  piv2 = vapply(x2, function(g) pivot %in% g$from, FUN.VALUE = F)

  # Sort pivot groups and put them first
  x1[piv1] = lapply(x1[piv1], sortGroup)
  x1 = x1[order(piv1, decreasing = T)]

  x2[piv2] = lapply(x2[piv2], sortGroup)
  x2 = x2[order(piv2, decreasing = T)]

  # Ensure locus 1 has max number of pivot groups (either 1 or 2)
  if(sum(piv1) < sum(piv2))
    x[] = list(locus1 = x2, locus2 = x1)
  else
    x[] = list(locus1 = x1, locus2 = x2)

  x
}


kinReduce = function(kin) {
  for(i in 1:2) for(j in seq_along(kin[[i]])) {
    g = kin[[i]][[j]]
    dups = duplicated(paste(g$from, g$to))
    if(any(dups))
      kin[[i]][[j]] = list(from = g$from[!dups], to = g$to[!dups])
  }
  kin
}

# Two locus kinship class
# Locus 1, replace `id` with par1 in group(s) gr1
# Locus 2, replace `id` with par2 in group(s) gr2
kinRepl = function(kin, id, loc1Rep, loc2Rep = NULL) {#par1, gr1 = 1, par2 = NULL, gr2 = 1) {

  kin$locus1 = kinRepl_1L(kin$locus1, id, loc1Rep$from1, loc1Rep$to1, loc1Rep$from2, loc1Rep$to2)

  if(!is.null(loc2Rep))
    kin$locus2 = kinRepl_1L(kin$locus2, id, loc2Rep$from1, loc2Rep$to1, loc2Rep$from2, loc2Rep$to2)

  validateKin2L(kin) # TODO: remove validate (coded to be unneccessary)
}

# Replacement at 1 locus
# G = list of kinship gorups ( = a pair of vectors (from, to) of same length)
# `id` is assumed to be present only in the first (one or two) group(s) of G
kinRepl_1L = function(G, id, from1, to1, from2 = NULL, to2 = NULL) {
  id = as.integer(id)
  from1 = as.integer(from1)
  to1 = as.integer(to1)
  from2 = as.integer(from2)
  to2 = as.integer(to2)
  if(length(id) != 1) stop2("Argument `id` must have length 1: ", id)
  if(length(from1) != length(to1)) stop2("Incompatible lengths of source vs indices: ", from1, " vs. ", to1)
  if(length(from2) != length(to2)) stop2("Incompatible lengths of source vs indices: ", from2, " vs. ", to2)
  if(anyNA(c(id, from1, to1, from2, to2))) stop2("NA's detected in recursion step: ", c(id, from1, to1, from2, to2))

  # Group 1: Replace id with `from1`, and corresponding entries in `to` with `to1`
  g = G[[1]]
  idx = g$from == id
  G[[1]] = list(from = c(from1, g$from[!idx]), to = c(to1, g$to[!idx]))

  # Group 2, if given
  if(length(from2)) {
    g = G[[2]]
    idx = g$from == id
    G[[2]] = list(from = c(from2, g$from[!idx]), to = c(to2, g$to[!idx]))
  }

  # Return modified G
  G
}

char2kinList = function(x) {
  # Example input: "5>9 = 7>10, 6>9, 8>10"
  # Output: list(from = c(5,7), to = c(9,10)), list(from = 6, to = 9), list(from = 8, to = 10)
  groups = strsplit(x, ",", fixed = T)[[1]]

  kinList = lapply(strsplit(groups, "=", fixed = T), function(g) {
    g = gsub("^[ ]*|[ ]*$", "", g)
    glist = strsplit(g, ">", fixed = T)

    from = sapply(glist, '[', 1)
    to = sapply(glist, '[', 2) # NA if no meiosis indicator

    list(from = as.integer(from), to = as.integer(to))
  })

  # Fix missing meiosis indicators
  USED.IDX = list()
  for(i in seq_along(kinList)) {
    from = kinList[[i]]$from; to = kinList[[i]]$to
    if(anyNA(to)) {
      for(j in which(is.na(to))) {
        a = as.character(from[j])
        prev = USED.IDX[[a]]
        val = if(is.null(prev)) 1 else prev + 1
        USED.IDX[[a]] = kinList[[i]]$to[j] = val
      }
    }
  }

  kinList
}

# Convert (single-locus) kinship group list to a string
kinList2char = function(x, labs = NULL) {
  if(is.null(labs)) {
    grvec = vapply(x, function(g) {
      paste(g$from, g$to, sep = ">", collapse = " = ")
    }, FUN.VALUE = "")
  }
  else {
    grvec = vapply(x, function(g) {
      from = labs[g$from]
      to = as.character(g$to)
      paste(from, to, sep=">", collapse=" = ")
    }, FUN.VALUE = "")
  }

  paste(grvec, collapse = ", ")
}


# Convert single-locus list of kinship groups to internal labels
kinList2internal = function(ped, kinList) {
  lapply(kinList, function(g) {
    g$from = internalID(ped, g$from)
    g$to = as.integer(g$to)
    return(g)

    #to = suppressWarnings(as.integer(g$to))
    #valid = is.na(to) | to > 0
    #to[valid] = internalID(ped, g$to[valid])
    #list(from = from, to = to)
  })
}

# Initialise memoisation used in twoLocusKinship and twoLocusGeneralisedKinship
initialiseTwoLocusMemo = function(ped, rho, recomb = NULL, chromType = "autosomal", counters = NULL) {
  if(chromType != "autosomal")
    stop2("Only `chromType = autosomal` is implemented at the moment")

  # Create memory storage
  mem = new.env()

  # Start timing
  st = Sys.time()

  mem$FIDX = ped$FIDX
  mem$MIDX = ped$MIDX
  mem$SEX = ped$SEX
  mem$rho = rho

  # Conditions on recombinant/non-recombinant gametes
  mem$recomb = mem$nonrecomb = rep(FALSE, pedsize(ped))
  if(!is.null(recomb)) {
    mem$recomb[internalID(ped, recomb$r)] = TRUE
    mem$nonrecomb[internalID(ped, recomb$nr)] = TRUE
  }

  # Logical matrix showing who has a common ancestor within the pedigree.
  # TODO: is this neccessary? (since we also include k1 below)
  mem$anc = hasCommonAncestor(ped)

  # Compute kinship matrix directly
  mem$k1 = switch(chromType, autosomal = kinship(ped), x = kinshipX(ped))

  # Storage for two-locus kinship values
  mem$k2 = list()

  # For quick look-up:
  FOU = founders(ped, internal = T)
  isFounder = rep(FALSE, pedsize(ped))
  isFounder[FOU] = TRUE
  mem$isFounder = isFounder

  ### Founder inbreeding
  # For linked loci only *completely* inbred founders are meaningful.
  finb = founderInbreeding(ped, ids = founders(ped), named = T)
  if(length(partials <- finb[finb > 0 & finb < 1]) > 0) {
    stop2("Partial founder inbreeding detected!",
          sprintf("\n  Individual '%s' (f = %g)", names(partials), partials),
          "\nTwo-locus coefficients are well-defined only when each founder is either outbred or completely inbred.")
  }

  # A logical vector for quick look-up e.g. isCompletelyInbred[a].
  isCompletelyInbred = rep(NA, pedsize(ped))
  isCompletelyInbred[FOU] = finb == 1
  mem$isCompletelyInbred = isCompletelyInbred

  # Counters
  for(cou in counters)
    assign(cou, 0, envir = mem)

  # Start time
  mem$st = st
  mem$initTime = Sys.time() - st

  mem
}
