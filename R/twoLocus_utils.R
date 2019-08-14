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
      # Check that meioses are compatible with pedigree
      target.in.ped = g$to > 0
      if(!is.null(ped) && length(target.in.ped) > 0) {
        labs = labels(ped)
        ch = g$to[target.in.ped]
        par = g$from[target.in.ped]
        par.ok = par == ped$FIDX[ch] | par == ped$MIDX[ch]
        if(!all(par.ok))
          printAndStop(x, toString(sprintf("%s is not a child of %s", labs[ch], labs[par])))
      }
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


sort_kin2L_OLD = function(x) {
  # Auxiliary function for sorting each group
  sortGroup = function(g) {
    ord = order(g$from, -g$to, decreasing = T)
    list(from = g$from[ord], to = g$to[ord])
  }

  x1 = lapply(x$locus1, sortGroup)
  x2 = lapply(x$locus2, sortGroup)

  # Sort on first 'from' of each group
  x1 = x1[order(vapply(x1, function(g) g$from[1], 1), decreasing = T)]
  x2 = x2[order(vapply(x2, function(g) g$from[1], 1), decreasing = T)]

  j1 = if (x1[[1]]$from[1] >= x2[[1]]$from[1]) 1 else 2
  j2 = 3-j1
  x[[j1]] = x1
  x[[j2]] = x2
  x
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
kinRepl = function(kin, id, par1, gr1 = 1, par2 = NULL, gr2 = 1) {
  y = kin
  id = as.integer(id)
  par1 = as.integer(par1)
  par2 = as.integer(par2)

  p1 = length(par1)
  p2 = length(par2)

  # Locus 1
  if(length(gr1) == 1 && p1 %in% 1:2) {
    G = y$locus1[[gr1]]
    mch = G$from == id
    id1 = rep_len(id, p1) # duplicate if par1 has length 2
    y$locus1[[gr1]] = list(from = c(par1, G$from[!mch]),
                           to = c(id1, G$to[!mch]))
  }
  else if(length(gr1) == 2 && p1 == 2) {
    for(i in 1:2) {
      G = y$locus1[[gr1[i]]]
      mch = G$from == id
      y$locus1[[gr1[i]]] = list(from = c(par1[i], G$from[!mch]),
                                to = c(id, G$to[!mch]))
    }
  }
  else stop2("Invalid input")

  if(p2 == 0)
    return(y)

  # Locus 2
  if(length(gr2) == 1 && p2 %in% 1:2) {
    G = y$locus2[[gr2]]
    mch = G$from == id
    id2 = rep_len(id, p2) # duplicate if par1 has length 2
    y$locus2[[gr2]] = list(from = c(par2, G$from[!mch]),
                           to = c(id2, G$to[!mch]))
  }
  else if(length(gr2) == 2 && p2 == 2) {
    for(i in 1:2) {
      G = y$locus2[[gr2[i]]]
      mch = G$from == id
      y$locus2[[gr2[i]]] = list(from = c(par2[i], G$from[!mch]),
                                to = c(id, G$to[!mch]))
    }
  }
  else stop2("Invalid input")

  validateKin2L(y)
}


char2kinList = function(x) {
  # Example input: "5>9 = 7>10, 6>9, 8>10"
  # Output: list(from = c(5,7), to = c(9,10)), list(from = 6, to = 9), list(from = 8, to = 10)
  groups = strsplit(x, ",", fixed = T)[[1]]

  kinList = lapply(strsplit(groups, "=", fixed = T), function(g) {
    g = gsub("^[ ]*|[ ]*$", "", g)
    glist = strsplit(g, ">", fixed = T)
    if(!all(lengths(glist) == 2))
      stop2("Something wrong with this kinship group: ", g)
    list(from = sapply(glist, '[', 1), to = sapply(glist, '[', 2))
  })

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
      valid = g$to > 0
      to[valid] = labs[g$to[valid]]
      paste(from, to, sep=">", collapse=" = ")
    }, FUN.VALUE = "")
  }

  paste(grvec, collapse = ", ")
}


char2kinList = function(x) {
  # Example input: "5>9 = 7>10, 6>9, 8>10"
  # Output: list(from = c(5,7), to = c(9,10)), list(from = 6, to = 9), list(from = 8, to = 10)
  groups = strsplit(x, ",", fixed = T)[[1]]

  kinList = lapply(strsplit(groups, "=", fixed = T), function(g) {
    g = gsub("^[ ]*|[ ]*$", "", g)
    glist = strsplit(g, ">", fixed = T)
    if(!all(lengths(glist) == 2))
      stop2("Something wrong with this kinship group: ", g)
    list(from = sapply(glist, '[', 1), to = sapply(glist, '[', 2))
  })

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
      valid = g$to > 0
      to[valid] = labs[g$to[valid]]
      paste(from, to, sep=">", collapse=" = ")
    }, FUN.VALUE = "")
  }

  paste(grvec, collapse = ", ")
}

# Convert single-locus list of kinship groups to internal labels
kinList2internal = function(ped, kinList) {
  lapply(kinList, function(g) {
    from = internalID(ped, g$from)
    to = suppressWarnings(as.integer(g$to))
    valid = is.na(to) | to > 0
    to[valid] = internalID(ped, g$to[valid])

    list(from = from, to = to)
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
