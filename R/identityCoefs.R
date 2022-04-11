#' Omnibus function for identity coefficients
#'
#' This function calculates the pairwise identity coefficients described by
#' Jacquard (1974). Unlike [condensedIdentity()] (and
#' `identity::identity.coefs()`), this function also computes the 15 *detailed*
#' identity coefficients. If `detailed = FALSE` the output is equivalent to
#' [condensedIdentity()], which is usually faster.
#'
#' Karigl (1981) gave the first recursive algorithm for the 9 condensed Jacqaurd
#' coefficients, based on certain *generalised kinship coefficients*. Karigl's
#' method is implemented in [condensedIdentity()].
#'
#' Weeks & Lange (1988) suggested a broader and more natural generalisation of
#' kinship coefficients, and a recursive algorithm for computing them. This is
#' implemented in `gKinship(..., method = "WL")`.
#'
#' Lange & Sinsheimer (1992) described an even further generalisation of kinship
#' coefficients, allowing a mix of deterministic and random sampling of alleles.
#' This is also implemented in `gKinship(..., method = "LS")`.
#'
#' In the same paper, Lange & Sinsheimer (1992) used the new generalisations to
#' give (i) an alternative algorithm for the 9 condensed identity coefficients,
#' and (ii) an algorithm for the 15 detailed coefficients. Both of these are
#' implemented in `identityCoefs()`.
#'
#' `method = "idcoefs"` calls the C program `IdCoefs` (version 2.1.1) by Mark
#' Abney (2009). This requires `IdCoefs` to installed on the computer (see link
#' under References) and available on the system search path.
#' `identity_idcoefs()` then writes the necessary files to disk and calls
#' `IdCoefs` via [system()].
#'
#' `method = "identity"` wraps `identity::identity.coefs()`, which is an R
#' interface for `IdCoefs`.
#'
#' Both the condensed and detailed coefficients are given in the orders used by
#' Jacquard (1974). The function `detailed2condensed` converts from detailed
#' coefficients (d1, ... d15) to condensed ones (D1, ..., D9) using the
#' following relations:
#'
#' * D1 = d1
#'
#' * D2 = d6
#'
#' * D3 = d2 + d3
#'
#' * D4 = d7
#'
#' * D5 = d4 + d5
#'
#' * D6 = d8
#'
#' * D7 = d9 + d12
#'
#' * D8 = d10 + d11 + d13 + d14
#'
#' * D9 = d15
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object.
#' @param ids A vector of two ID labels.
#' @param detailed A logical. If FALSE (default), the 9 condensed coefficients
#'   are computed; otherwise the 15 detailed identity coefficients.
#' @param Xchrom A logical, by default FALSE.
#' @param self A logical indicating if self-relationships (e.g., between a
#'   pedigree member and itself) should be included. FALSE by default.
#' @param method Either "auto", "K", "WL", "LS", "GC", "identity" or "idcoefs".
#' @param simplify Simplify the output (to a numeric of length 9) if `ids` has
#'   length 2. Default: TRUE.
#' @param verbose A logical.
#' @param ... Further arguments.
#' @param d Either a numeric vector of length 15, or a data frame with 17
#'   columns.
#'
#' @return A data frame with L + 2 columns, where L is either 9 (default) or 15
#'   (if `detailed = TRUE`).
#'
#'   If `simplify = TRUE` and `length(ids) = 2`: A numeric vector of length `L`.
#'
#' @seealso [condensedIdentity()], [gKinship()]
#'
#' @references
#'
#' * Jacquard, A. (1974). The Genetic Structure of Populations. Springer.
#'
#' * Karigl, G. (1981). A recursive algorithm for the calculation of identity
#' coefficients. Ann. Hum. Genet.
#'
#' * Weeks, D.E. & Lange, K. (1988). The affected-pedigree-member method of
#' linkage analysis. Am. J. Hum. Genet
#'
#' * Lange, K. & Sinsheimer, J.s. (1992). Calculation of genetic identity
#' coefficients. Ann. Hum. Genet.
#'
#' * Abney, Mark (2009). A graphical algorithm for fast computation of identity
#' coefficients and generalized kinship coefficients. Bioinformatics, 25,
#' 1561-1563. <https://home.uchicago.edu/~abney/abney_web/Software.html>

#' @examples
#' x = fullSibMating(1)
#' ids = 1:6
#'
#' ### Condensed coefficients
#' j1 = identityCoefs(x, ids, method = "K")
#' j2 = identityCoefs(x, ids, method = "WL")
#' j3 = identityCoefs(x, ids, method = "LS")
#' j4 = identityCoefs(x, ids, method = "GC")
#' j5 = condensedIdentity(x, ids) # legacy version
#' # j6 = identityCoefs(x, ids, method = "idcoefs") # requires IdCoefs installed
#'
#' stopifnot(all.equal(j1,j2), all.equal(j1,j3), all.equal(j1,j4), all.equal(j1,j5))
#'
#' ### Detailed coefficients
#' jdet1 = identityCoefs(x, ids, detailed = TRUE, method = "LS")
#' jdet2 = identityCoefs(x, ids, detailed = TRUE, method = "GC")
#' stopifnot(all.equal(jdet1,jdet2))
#'
#' ### X chromosomal coefficients
#' jx1 = identityCoefs(x, ids, method = "K", Xchrom = TRUE)
#' jx2 = identityCoefs(x, ids, method = "GC", Xchrom = TRUE)
#' jx3 = condensedIdentityX(x, ids)  # legacy
#' stopifnot(all.equal(jx1,jx2), all.equal(jx1,jx3))
#'
#' @export
identityCoefs = function(x, ids, detailed = FALSE, Xchrom = FALSE, self = FALSE, simplify = TRUE,
                    method = c("auto", "K", "WL", "LS", "GC", "identity", "idcoefs"),
                    verbose = FALSE, ...) {

  method = match.arg(method)
  if(method == "auto")
    method = chooseIdentityMethod(x, ids, detailed, Xchrom)

  res = switch(method,
    K = {
      if(detailed)
        stop2("Method 'K' does not support detailed coefficients. Try e.g., 'LS' or 'GC'.")
      identity_Karigl(x, ids, Xchrom = Xchrom, self = self, verbose = verbose, ...)
      },
    WL = {
      if(detailed)
        stop2("Method 'WL' does not support detailed coefficients. Try e.g., 'K', 'LS' or 'GC'.")
      identity_WL(x, ids, Xchrom = Xchrom, self = self, verbose = verbose)
    },
    LS = {
      if(Xchrom)
        stop2("Method 'LS' does not support X-chromosomal coefficients. Try e.g., 'K'.")
      identity_LS(x, ids, detailed = detailed, self = self, verbose = verbose)
      },
    GC = {
      identity_GC(x, ids, detailed = detailed, Xchrom = Xchrom, self = self, verbose = verbose)
      },
    identity = {
      if(Xchrom)
        stop2("The 'identity' package does not support X-chromosomal coefficients.")
      if(detailed)
        stop2("The 'identity' package does not support detailed coefficients.")
      identity_identity(x, ids, self = self, verbose = verbose)
    },
    idcoefs = {
      if(Xchrom)
        stop2("The 'IdCoefs' program does not support X-chromosomal coefficients.")
      if(detailed)
        stop2("The 'IdCoefs' program does not support detailed coefficients.")
      identity_idcoefs(x, ids, self = self, verbose = verbose, ...)
    }
  )

  # For single pair: By default simplify to unnamed vector
  if(nrow(res) == 1 && simplify)
    return(as.numeric(res[-(1:2)]))

  res
}


prepPed = function(x, addpar = NULL, Xchrom = FALSE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  if(!is.null(addpar))
    x = addFounderParents(x, ids = addpar, Xchrom = Xchrom)

  # Ensure all founders are put first
  x = foundersFirst(x)

  x
}

prepIds = function(x, ids, self = FALSE, int = FALSE) {

  # Special case (for convenience) regardless of `self`
  if(length(ids) == 2 && ids[1] == ids[2]) {
    if(int) ids = internalID(x, ids)
    return(list(ids))
  }

  ids = as.character(unique.default(ids))
  idsInt = internalID(x, ids)

  if(int)
    ids = idsInt

  # Always sort on internal
  ids =  .mysort(ids, by = idsInt)

  # List of all pairs
  .idPairs(ids, self = self, as = if(int) "integer" else "character")
}


# Setup memoisation
memoIdentity = function(x, Xchrom = FALSE, method = NULL, counters = NULL, maxId = pedsize(x), sparse = 30, verbose = FALSE) {
  mem = new.env()
  mem$st = Sys.time()
  mem$method = method
  mem$FIDX = x$FIDX
  mem$MIDX = x$MIDX
  mem$SEX = x$SEX

  # Kinship matrix and inbreeding coefficients
  mem$KIN = kinship(x, Xchrom = Xchrom)
  mem$INB = 2*diag(mem$KIN) - 1
  if(Xchrom)
    mem$INB[getSex(x) == 1L] = 1 # on X, males are always 1

  # Related pairs: REL[i,j] = TRUE iff are related
  mem$REL = mem$KIN > 0

  # Storage for generalised gamete probabilities
  mem$MEM = list()

  # For quick look-up of founders
  mem$isFounder = labels(x) %in% founders(x)

  # Counters
  counters = switch(method,
    K = c("i2", "i3", "i4", "i22", "i2r", "i3r", "i4r", "i22r"),
    WL = c("i", "B0", "B1", "B2", "B3", "ilook", "irec"),
    LS = c("i", "B0", "B1", "B2", "B3", "D1", "D2", "ilook", "irec"),
    GC = c("i", "ilook", "irec", "i22")
  )

  for(cou in counters)
    assign(cou, value = 0L, envir = mem)

  ### Special stuff for some methods
  if(method == "LS") {
    # Second storage for deterministic patterns
    mem$MEM2 = list()
  }
  if(method == "K") {
    mem$ANC = mem$REL

    # When to use sparse arrays. Default: pedsize > 30
    if(length(sparse) != 1 || is.character(sparse))
      stop2("Argument `sparse` must be a number or a logical: ", sparse)
    sparse = if(is.numeric(sparse)) maxId > sparse else isTRUE(sparse)

    if(sparse) {
      if(verbose) message("Using sparse lookup tables")
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

    for(i in founders(x, internal = TRUE)) {
      if(i > maxId) break # otherwise out of range!
      fi = mem$INB[i]
      KIN3[i, i, i] = (1 + 3*fi)/4
      KIN4[i, i, i, i] = (1 + 7*fi)/8
      KIN22[i, i, i, i] = (1 + 3*fi)/4
    }

    mem$KIN3 = KIN3
    mem$KIN4 = KIN4
    mem$KIN22 = KIN22
  }
  mem
}

chooseIdentityMethod = function(x, ids, detailed, Xchrom) {
  if(!detailed)
    return("K")

  # Detailed on X: only GC
  if(Xchrom)
    return("GC")

  # Detailed autosomal: LS is fastest, but requires adding parents to founders.
  # If inbred founders involved: GC
  hasIF = hasInbredFounders(x) && any(inbreeding(x)[intersect(ids, founders(x))] > 0)
  if(hasIF)
    return("GC")

  return("LS")
}

printMemInfo = function(mem) {
  tot_time = format(Sys.time()-mem$st, digits = 3)

  msg = switch(mem$method,
    LS = {
      b = c(mem$B0, mem$B1, mem$B2, mem$B3)
      d = c(mem$D1 %||% 0, mem$D2 %||% 0)
      glue::glue("
                 Function calls: {mem$i}
                   Boundary:     {sum(b)} (B0={b[1]}, B1={b[2]}, B2={b[3]}, B3={b[4]})
                   Determ:       {sum(d)} (D1={d[1]}, D2={d[2]})
                   Lookups:      {mem$ilook}
                   Recursions:   {mem$irec}
                 Total time: {tot_time}")
    },
    WL = {
      b = c(mem$B0, mem$B1, mem$B2, mem$B3)
      glue::glue("
                 Function calls: {mem$i}
                   Boundary:     {sum(b)} (B0={b[1]}, B1={b[2]}, B2={b[3]}, B3={b[4]})
                   Lookups:      {mem$ilook}
                   Recursions:   {mem$irec}
                 Total time: {tot_time}")
    },
    GC = {
      glue::glue("
                Function calls:
                   psi2: {mem$i[2]} ({mem$ilook[2]} lookups, {mem$irec[2]} recursions
                   psi3: {mem$i[3]} ({mem$ilook[3]} lookups, {mem$irec[3]} recursions
                   psi4: {mem$i[4]} ({mem$ilook[4]} lookups, {mem$irec[4]} recursions
                  psi22: {mem$i22} ({mem$ilook22} lookups, {mem$irec22} recursions
                Total time used: {tot_time}")
    },
    K = {
      glue::glue("
                   Function calls:
                      phi2: {mem$i2} (recursions: {mem$i2r}
                      phi3: {mem$i3} (recursions: {mem$i3r})
                      phi4: {mem$i4} (recursions: {mem$i4r})
                     phi22: {mem$i22} (recursions: {mem$i22r})
                   Total time used: {tot_time}")
    })
  message(msg)
}



#' @rdname identityCoefs
#' @export
detailed2condensed = function(d) {
  ok = (is.data.frame(d) && ncol(d) == 17) || (is.numeric(d) && length(d) == 15)
  if(!ok)
    stop2("Input must be either a vector of length 15, or a data frame with 17 columns")

  if(is.data.frame(d)) {

    idcols = d[, 1:2]
    d = d[,-(1:2)] # coefficient columns

    return(cbind(idcols,
      D1 = d[[1]],
      D2 = d[[6]],
      D3 = d[[2]] + d[[3]],
      D4 = d[[7]],
      D5 = d[[4]] + d[[5]],
      D6 = d[[8]],
      D7 = d[[9]] + d[[12]],
      D8 = d[[10]] + d[[11]] + d[[13]] + d[[14]],
      D9 = d[[15]]))
  }

  # Else: Numeric vector
  c(d[[1]],
    d[[6]],
    d[[2]] + d[[3]],
    d[[7]],
    d[[4]] + d[[5]],
    d[[8]],
    d[[9]] + d[[12]],
    d[[10]] + d[[11]] + d[[13]] + d[[14]],
    d[[15]]
  )
}

# For X-chromosomal coefficients, put NAs in ill-defined entries
Xmask = function(x, df) {
  if(!is.data.frame(df) || !ncol(df) %in% c(11,17))
    stop2("The `df` argument must be a data frame with either 11 or 17 columns")
  detailed = ncol(df) == 17

  sex1 = getSex(x, df[,1])
  sex2 = getSex(x, df[,2])
  if(detailed) {
    df[sex1 == 1 & sex2 == 1, 2 + c(2,3,4,5,7:15)] = NA
    df[sex1 == 1 & sex2 == 2, 2 + c(4,5,8,9:15)] = NA
    df[sex1 == 2 & sex2 == 1, 2 + c(2,3,7,9:15)] = NA
  }
  else {
    df[sex1 == 1 & sex2 == 1, 2 + 3:9] = NA
    df[sex1 == 1 & sex2 == 2, 2 + 5:9] = NA
    df[sex1 == 2 & sex2 == 1, 2 + c(3:4,7:9)] = NA
  }
  df
}


debugReturn = function(res, debug = FALSE, indent = 0, comment = NULL) {
  if(debug)
    message(strrep(" ", indent), res, comment)
  res
}
