initialiseMemo = function(ped, ids, sparse = 20, chromType = "autosomal", verbose = FALSE) {

  chromType = match.arg(tolower(chromType), c("autosomal", "x"))

  # Create memory storage
  mem = new.env()

  # Start timing
  st = Sys.time()

  FIDX = ped$FIDX
  MIDX = ped$MIDX
  SEX = ped$SEX

  # Logical matrix showing who has a common ancestor within the pedigree.
  anc = hasCommonAncestor(ped)

  # Compute kinship matrix directly
  KIN2 = switch(chromType, autosomal = kinship(ped), x = kinshipX(ped))

  maxId = max(ids)
  if(is.na(sparse)) sparse = 20
  if(is.numeric(sparse)) sparse = maxId > sparse

  if(sparse) {
    if(verbose) cat("Using sparse lookup tables\n")
    sparsarr0 = slam::simple_sparse_zero_array
    KIN3 = sparsarr0(dim = rep(maxId, 3), mode = "double")
    KIN4 = sparsarr0(dim = rep(maxId, 4), mode = "double")
    KIN22 = sparsarr0(dim = rep(maxId, 4), mode = "double")
  }
  else {
    KIN3 = array(NA_real_, dim = rep(maxId, 3))
    KIN4 = array(NA_real_, dim = rep(maxId, 4))
    KIN22 = array(NA_real_, dim = rep(maxId, 4))
  }

  # For quick look-up:
  FOU = founders(ped, internal = T)
  isFounder = rep(FALSE, pedsize(ped))
  isFounder[FOU] = TRUE

  # Founder inbreeding
  # A vector of length pedsize(ped), with inb.coeffs at all founder idx,
  # and NA entries everywhere else. Enables quick look-up e.g. founderInb[a].
  founderInb = rep(NA_real_, pedsize(ped))
  founderInb[FOU] = founderInbreeding(ped, ids = founders(ped), chromType = chromType)

  for(i in FOU) {
    if(i > maxId) break # otherwise out of range!
    fi = founderInb[i]
    KIN3[i, i, i] = (1 + 3*fi)/4
    KIN4[i, i, i, i] = (1 + 7*fi)/8
    KIN22[i, i, i, i] = (1 + 3*fi)/4
  }

  mem$ANC = anc
  mem$FIDX = FIDX
  mem$MIDX = MIDX
  mem$SEX = SEX

  mem$isFounder = isFounder
  mem$founderInb = founderInb

  mem$KIN2 = KIN2
  mem$KIN3 = KIN3
  mem$KIN4 = KIN4
  mem$KIN22 = KIN22

  # Counters
  mem$i2 = mem$i3 = mem$i4 = mem$i22 = 0
  mem$i2r = mem$i3r = mem$i4r = mem$i22r = 0

  # Start time
  mem$st = st
  mem$initTime = Sys.time() - st

  mem
}

printCounts = function(mem) {
  init_time = format(mem$initTime, digits = 4)
  tot_time = format(Sys.time()-mem$st, digits = 4)

  msg = glue::glue("
                   Function calls:
                   phi2  = {mem$i2}
                   phi3  = {mem$i3} (recursions: {mem$i3r})
                   phi4  = {mem$i4} (recursions: {mem$i4r})
                   phi22 = {mem$i22} (recursions: {mem$i22r})
                   Total time used: {tot_time}")
  print(msg)
}
