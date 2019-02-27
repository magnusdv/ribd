#' Identity coefficients on X
#'
#' Computes the X chromosomal condensed identity coefficients of a pairwise
#' relationship.
#'
#' The implementation is inspired by Karigl's recursive algorithm (1981) for the
#' autosomal case, modified to account for X-linked inheritance.
#'
#'
#' @param x A pedigree object
#' @param ids A numeric of length 2 containing ID labels of two pedigree members
#' @param verbose A Logical
#' @param checkAnswer If TRUE, the result is checked against the output of the
#'   XIBD package, if this is installed. (Ignored if any of the founders are inbred.)
#' @param sparse A positive integer, indicating the pedigree size limit for
#'   using sparse arrays (as implemented by the
#'   [slam](https://CRAN.R-project.org/package=slam) package) instead of
#'   ordinary arrays.
#'
#' @return A vector of length 9. Unless both of `ids` are female, some entries
#'   are NA.
#'
#' @export
#'
#' @examples
#' library(pedtools)
#'
#' x = fullSibMating(1)
#' x_sisters = swapSex(x, 5)
#' x_brothers = swapSex(x, 6)
#'
#' condensedIdentityX(x, ids = 5:6)
#' condensedIdentityX(x_sisters, ids = 5:6)
#' condensedIdentityX(x_brothers, ids = 5:6)
#'
condensedIdentityX = function(x, ids, verbose=FALSE, checkAnswer=verbose, sparse=Inf) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must be a vector of length 2")
  if(any(founderInbreeding(x) > 0)) stop2("Inbred founders are not yet implemented for this function")
  if(sparse < Inf) stop2("Sparse arrays are not yet implemented for this function.")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  FIDX = x$FIDX
  MIDX = x$MIDX
  SEX = x$SEX

  phi2 = function(a, b) {
    i2 <<- i2+1
    if(a*b == 0) return(0)

    # Sort: a>=b
    if(a < b) {tmp=a; a=b; b=tmp}

    if(SEX[a] == 1) {
      if(a == b) 1
      else phi2(MIDX[a], b)
    }
    else {
      if(a==b) (1 + phi2(FIDX[a], MIDX[a]))/2
      else (phi2(FIDX[a], b) + phi2(MIDX[a], b))/2
    }
  }

  phi3 = function(a, b, c) {
    i3 <<- i3+1
    if(a*b*c == 0) return(0)

    # Sort: a>=b>=c
    if(a < b) {tmp=a; a=b; b=tmp}
    if(b < c) {tmp=b; b=c; c=tmp}
    if(a < b) {tmp=a; a=b; b=tmp}

    if(SEX[a] == 1) {
      if(a == b && a == c) 1
      else if(a == b) phi2(MIDX[a], c)
      else phi3(MIDX[a], b, c)
    }
    else {
      if(a == b && a == c)
        (1 + 3*phi2(FIDX[a], MIDX[a]))/4
      else if(a == b)
        (phi2(a, c) + phi3(FIDX[a], MIDX[a], c))/2
      else
        (phi3(FIDX[a], b, c) + phi3(MIDX[a], b, c))/2
    }
  }

  phi4 = function(a, b, c, d) {
    i4 <<- i4+1
    if(a*b*c*d == 0) return(0)

    # Sort: a >= b >= c >= d
    s = sort.int(c(a,b,c,d),method="quick")
    a = s[4]; b = s[3]; c = s[2]; d = s[1]

    if(SEX[a] == 1) {
      if(a == b && a == c && a == d)
        1
      else if(a == b && a == c)
        phi2(a, d)
      else if(a == b)
        phi3(MIDX[a], c, d)
      else
        phi4(MIDX[a], b, c, d)
    }
    else {
      if(a == b && a == c && a == d)
        (1 + 7*phi2(FIDX[a], MIDX[a]))/8
      else if(a == b && a == c)
        (phi2(a, d) + 3*phi3(FIDX[a], MIDX[a], d))/4
      else if(a == b)
        (phi3(a, c, d) + phi4(FIDX[a], MIDX[a], c, d))/2
      else
        (phi4(FIDX[a], b, c, d) + phi4(MIDX[a], b, c, d))/2
    }
  }

  phi22 = function(a, b, c, d) {
    i22 <<- i22+1
    if(a*b * c*d == 0) return(0) # Karigl-81: a*b == 0 AND c*d == 0. Misprint??

    # Sort: a >= b,c,d; c >= d; if(a==c) then b>=d
    s = c(min(a,b), max(a,b), min(c,d), max(c,d)) # d,c,b,a
    if(s[4] < s[2] || (s[4] == s[2] && s[3] < s[1]))
      s[] = s[c(3,4,1,2)]
    a = s[4]; b = s[3]; c = s[2]; d = s[1]

    if(SEX[a] == 1) {
      if(a == b && a == c && a == d)
        1
      else if(a == b && a == c)
        phi2(MIDX[a], d)
      else if(a == b)
        phi2(c, d)
      else if(a == c)
        phi3(MIDX[a], b, d)
      else
        phi22(MIDX[a], b, c, d)
    }
    else {
      if(a == b && a == c && a == d)
        (1 + 3*phi2(FIDX[a], MIDX[a]))/4
      else if(a == b && a == c)
        (phi2(a, d) + phi3(FIDX[a], MIDX[a], d))/2
      else if(a == b)
        (phi2(c, d) + phi22(FIDX[a], MIDX[a], c, d))/2
      else if(a == c)
        (2*phi3(a, b, d) + phi22(FIDX[a], b, MIDX[a], d) + phi22(MIDX[a], b, FIDX[a], d))/4
      else
        (phi22(FIDX[a], b, c, d) + phi22(MIDX[a], b, c, d))/2
    }
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
    16,0,4,0,4,0,4,1,0), byrow=TRUE, ncol=9)

  st = Sys.time()

  # Initialize counters
  i2 <- i3 <- i4 <- i22 <- 0

  id1 = ids_int[1]; id2 = ids_int[2]
  RHS = c(1,
          2*phi2(id1,id1),
          2*phi2(id2,id2),
          4*phi2(id1,id2),
          8*phi3(id1,id1,id2),
          8*phi3(id1,id2,id2),
          16*phi4(id1,id1,id2,id2),
          4*phi22(id1,id1,id2,id2),
          16*phi22(id1,id2,id1,id2))

  j = solve(M9, RHS)
  if(SEX[id1] == 1 && SEX[id2] == 1)
    j[3:9] = NA
  if(SEX[id1] == 1 && SEX[id2] == 2)
    j[5:9] = NA
  if(SEX[id1] == 2 && SEX[id2] == 1)
    j[c(3:4,7:9)] = NA


  if(verbose) {
    secs = sprintf("%.1f", Sys.time()-st)
    msg = glue::glue("
    Function calls:
      phi2  = {i2}
      phi3  = {i3}
      phi4  = {i4}
      phi22 = {i22}
    Time used: {secs} seconds")
    print(msg)
  }

  if(checkAnswer) compare_with_XIBD(x, ids, j)

  j
}

compare_with_XIBD = function(x, ids, j) {
  cat("Comparison with `XIBD` package: ")

  if(any(founderInbreeding(x) > 0)) {
    message("skipped. (Pedigree has inbred founders.)")
    return()
  }

  jj = xibd(x, ids)
  if(is.null(jj)) {
    message("Install `XIBD` or use `checkAnswer = FALSE` to avoid this message.")
    return()
  }

  if(identical(j, jj))
    message("OK!")
  else if(isTRUE(all.equal(j,jj)))
    message("all.equal() OK, but not identical()")
  else {
    message("*** MISMATCH! ***")
    cat("IDS:", ids, "\n")
    print(rbind(`XIBD:` = jj, `ribd:` = j))
  }
}

#' @importFrom utils capture.output
xibd = function(x, ids) {
  if(!requireNamespace("XIBD", quietly = TRUE)){
    message("Package `XIBD` is not installed.")
    return()
  }

  sex1 = getSex(x, ids[1])
  sex2 = getSex(x, ids[2])
  ids_int = internalID(x, ids)

  pedm = as.matrix(x)[, 1:4]
  colnames(pedm) = c("iid","pid", "mid","sex")

  capture.output(
    deltas <- lapply(1:9, function(i) XIBD:::delta(pedm, i, ids_int[1], ids_int[2]))
  )

  # If any of ids are male, the `deltas` list contains NULL entries
  jx = rep(NA_real_, 9)

  if(sex1 == 1 && sex2 == 1)
    jx[1:2] = unlist(deltas)
  else if (sex1 == 1 && sex2 == 2)
    jx[1:4] = unlist(deltas)
  else if (sex1 == 2 && sex2 == 1)
    jx[c(1,2,5,6)] = unlist(deltas)
  else
    jx[] = unlist(deltas)

  jx
}


