# Initialise memoisation used in twoLocusKinship and twoLocusGeneralisedKinship
initialiseTwoLocusMemo = function(ped, rho, recomb = NULL, chromType = "autosomal", counters = NULL) {
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
  for(cou in counters)
    assign(cou, 0, envir = mem)

  #mem$i = mem$ilook = mem$irec = 0
  #mem$eq7 = mem$eq8 = mem$eq9a = mem$eq9b = mem$eq10 = mem$eq11a = mem$eq11b = 0

  # Start time
  mem$st = st
  mem$initTime = Sys.time() - st

  mem
}
