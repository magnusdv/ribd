#' Identity coefficients
#'
#' Computes the 9 condensed identity coefficients of a pairwise relationship.
#' Founders of the pedigree may be inbred; use [pedtools::founder_inbreeding()]
#' to set this up.
#'
#' The implementation is a modified version of Karigl's recursive algorithm
#' (1981).
#'
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#' @param ids A numeric of length 2 containing ID labels of two pedigree
#'   members.
#' @param verbose A Logical
#' @param checkAnswer If TRUE, and the `identity` package is installed, the
#'   result is checked against the output of [identity::identity.coefs()].
#'   (Ignored if any of the founders are inbred.)
#' @param sparse A positive integer, indicating the pedigree size limit for
#'   using sparse arrays (as implemented by the
#'   [slam](https://CRAN.R-project.org/package=slam) package) instead of
#'   ordinary arrays.
#'
#' @return A vector of length 9, containing the condensed identity coefficients.
#' @references G. Karigl (1981). _A recursive algorithm for the calculation of
#'   identity coefficients_ Annals of Human Genetics, vol. 45.
#'
#' @seealso [pedtools::ped()], [pedtools::founder_inbreeding()]
#' @export
#'
#' @examples
#' library(pedtools)
#'
#' x = fullSibMating(2)
#' j1 = ibd_identity(x, ids = 5:6)
#'
#' stopifnot(all.equal(j1, c(2, 1,4, 1, 4, 1, 7, 10, 2)/32))
#'
#' # Recalculate the coefficients when individual 1 is 100% inbred
#' founder_inbreeding(x, 1) = 1
#' ibd_identity(x, ids = 5:6)
#'
#'
ibd_identity = function(x, ids, verbose=TRUE, checkAnswer=verbose, sparse=50) {
  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  FID = x$FID
  MID = x$MID
  FOU = founders(x, internal = T)

  # For quick look-up:
  is_founder = rep(FALSE, pedsize(x))
  is_founder[FOU] = TRUE

  # Logical matrix showing who has a common ancestor within the pedigree.
  anc = hasCA(x)

  # Compute kinship matrix directly
  KIN2 = ibd_kinship(x)

  mx_id = max(ids_int)
  use_sparse = mx_id > sparse

  if(use_sparse) {
    if(verbose) cat("Using sparse lookup tables\n")
    sparsarr0 = slam::simple_sparse_zero_array
    KIN3 = sparsarr0(dim = rep(mx_id, 3), mode = "integer")
    KIN4 = sparsarr0(dim = rep(mx_id, 4), mode = "integer")
    KIN22 = sparsarr0(dim = rep(mx_id, 4), mode = "integer")
  }
  else {
    KIN3 = array(NA, dim = rep(mx_id, 3))
    KIN4 = array(NA, dim = rep(mx_id, 4))
    KIN22 = array(NA, dim = rep(mx_id, 4))
  }

  # Founder inbreeding
  # A vector of length pedsize(x), with inb.coeffs at all founder idx,
  # and NA entries everywhere else. Enables quick look-up e.g. FOU_INB[a].
  FOU_INB = rep(NA_real_, pedsize(x))
  FOU_INB[FOU] = founder_inbreeding(x, ids=founders(x))

  for(i in FOU) {
    if(i > mx_id) break # otherwise out of range!
    fi = FOU_INB[i]
    KIN3[i, i, i] = (1 + 3*fi)/4
    KIN4[i, i, i, i] = (1 + 7*fi)/8
    KIN22[i, i, i, i] = (1 + 3*fi)/4
  }

  phi2 = function(a, b) {
    i2 <<- i2+1
    if(a*b == 0) return(0)
    KIN2[[a, b]]
  }

  phi3 = function(a, b, c) {
    i3 <<- i3+1
    if(a*b*c == 0) return(0)
    if(!all(anc[a,b],anc[b,c],anc[a,c])) return(0)

    # Sort: a>=b>=c
    if(a < b) {tmp=a; a=b; b=tmp}
    if(b < c) {tmp=b; b=c; c=tmp}
    if(a < b) {tmp=a; a=b; b=tmp}

    # Recursion function (called only if necessary)
    # Assumes a,b,c sorted
    phi3_recurse = function(a, b, c) {
      i3r <<- i3r + 1
      if(a == b && a == c)
        return((1 + 3*phi2(FID[a], MID[a]))/4)

      if(a == b)
        return((phi2(a, c) + phi3(FID[a], MID[a], c))/2)

      return((phi3(FID[a], b, c) + phi3(MID[a], b, c))/2)
    }

    # Lookup in array; compute if necessary.
    if(is.na(KIN3[[a,b,c]]))
      KIN3[[a,b,c]] <<- phi3_recurse(a, b, c)

    KIN3[[a,b,c]]
  }

  phi4 = function(a, b, c, d) {
    i4 <<- i4+1
    if(a*b*c*d == 0) return(0)
    if(!all(anc[a,b],anc[b,c],anc[c,d],anc[b,d],anc[a,c],anc[a,d])) return(0)

    # Sort: a >= b >= c >= d
    if(a < b) {tmp=a; a=b; b=tmp}
    if(a < c) {tmp=a; a=c; c=tmp}
    if(a < d) {tmp=a; a=d; d=tmp}
    if(b < c) {tmp=b; b=c; c=tmp}
    if(b < d) {tmp=b; b=d; d=tmp}
    if(c < d) {tmp=c; c=d; d=tmp}

    # Recursion function (called only if necessary)
    # Assumes a,b,c,d sorted
    phi4_recurse = function(a, b, c, d) {
      if(a == b && a == c && a == d)
        return((1 + 7*phi2(FID[a], MID[a]))/8)

      if(a == b && a == c)
        return((phi2(a, d) + 3*phi3(FID[a], MID[a], d))/4)

      if(a == b)
        return((phi3(a, c, d) + phi4(FID[a], MID[a], c, d))/2)

      return((phi4(FID[a], b, c, d) + phi4(MID[a], b, c, d))/2)
    }

    # Lookup in array; compute if necessary.
    if(is.na(KIN4[[a,b,c,d]]))
      KIN4[[a,b,c,d]] <<- phi4_recurse(a, b, c, d)

    KIN4[[a,b,c,d]]
  }

  phi22 = function(a, b, c, d) {
    i22 <<- i22+1
    if(a*b*c*d == 0) return(0)
    if(!(anc[a,b] && anc[c,d])) return(0)

    # Sort: a >= b,c,d; c >= d; if(a==c) then b>=d
    s = c(min(a,b), max(a,b), min(c,d), max(c,d)) # d,c,b,a
    if(s[4] < s[2] || (s[4] == s[2] && s[3] < s[1]))
      s[] = s[c(3,4,1,2)]
    a = s[4]; b = s[3]; c = s[2]; d = s[1]

    # Recursion function (called only if necessary)
    # Assumes a,b,c sorted
    phi22_recurse = function(a, b, c, d) {
      i22r <<- i22r+1
      if(a == b && a == c && a == d)
        return((1 + 3*phi2(FID[a], MID[a]))/4)

      if(a == b && a == c)
        return((phi2(a, d) + phi3(FID[a], MID[a], d))/2)

      if(a == b) { #NB modification to allow inbred founders!
        if(is_founder[a])
          return(0.5*phi2(c, d)*(1 + FOU_INB[a]))

        return((phi2(c, d) + phi22(FID[a], MID[a], c, d))/2)
      }

      if(a == c)
        return((2*phi3(a, b, d) + phi22(FID[a], b, MID[a], d) +
                  phi22(MID[a], b, FID[a], d))/4)

      return((phi22(FID[a], b, c, d) + phi22(MID[a], b, c, d))/2)
    }

    # Lookup in array; compute if necessary.
    if(is.na(KIN22[[a,b,c,d]]))
      KIN22[[a,b,c,d]] <<- phi22_recurse(a, b, c, d)

    KIN22[[a,b,c,d]]
  }

  M9 = matrix(c(
    1,1,1,1,1,1,1,1,1,
    2,2,2,2,1,1,1,1,1,
    2,2,1,1,2,2,1,1,1,
    4,0,2,0,2,0,2,1,0,
    8,0,4,0,2,0,2,1,0,
    8,0,2,0,4,0,2,1,0,
    16,0,4,0,4,0,2,1,0,
    4,4,2,2,2,2,1,1,1,
    16,0,4,0,4,0,4,1,0), byrow=T, ncol=9)

  st = Sys.time()

  # Initialize counters
  i2 <- i3 <- i3r <- i4 <- i22 <-i22r <- 0

  id1 = ids_int[1]; id2 = ids_int[2]
  RHS = c(
    1,
    2*phi2(id1,id1),
    2*phi2(id2,id2),
    4*phi2(id1,id2),
    8*phi3(id1,id1,id2),
    8*phi3(id1,id2,id2),
    16*phi4(id1,id1,id2,id2),
    4*phi22(id1,id1,id2,id2),
    16*phi22(id1,id2,id1,id2))

  #print(RHS)
  j = solve(M9, RHS)

  if(verbose) {
    secs = sprintf("%.1f", Sys.time()-st)
    msg = glue::glue("
    Function calls:
      phi2  = {i2}
      phi3  = {i3} (recurse: {i3r})
      phi4  = {i4}
      phi22 = {i22} (recurse: {i22r})
    Time used: {secs} seconds")
    print(msg)
  }

  if(checkAnswer) compare_with_identity(x, ids, j)

  j
}

compare_with_identity = function(x, ids, j) {
  cat("Comparison with `identity` package: ")

  if(any(founder_inbreeding(x) > 0)) {
    message("skipped. (Pedigree has inbred founders.)")
    return()
  }

  if(!requireNamespace("identity", quietly = TRUE)) {
    message("skipped. Package `identity` is not installed.")
    return()
  }

  jj = jacquard(x, ids)
  if(isTRUE(all.equal(j, jj)))
    message("OK!")
  else {
    message("*** MISMATCH! ***")
    cat("IDS:", ids, "\n")
    print(rbind(`identity:` = jj, `ribd:` = j))
  }
}

# Define double bracket extract/replace operators to accommodate the impossibility of zero-values (-1 used instead)
`[[.simple_sparse_array` = function(x, ...) {
  val = x[...]$v
  if(length(val) == 0) return(NA)
  if(val < 0) return(0)
  val
}

`[[<-.simple_sparse_array` = function(x, ..., value) {
  if(value == 0) value = -1
  x[...] <- value
  x
}
