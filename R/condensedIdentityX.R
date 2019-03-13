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
condensedIdentityX = function(x, ids, sparse = NA, verbose = FALSE, checkAnswer = verbose) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # Enforce parents to precede their children
  if(!has_parents_before_children(x))
    x = parents_before_children(x)

  ids_int = internalID(x, ids)

  # Setup memoisation
  mem = initialiseMemo(x, ids_int, sparse = sparse, chromType = "x", verbose = verbose)

  M9 = matrix(c(
    1,1,1,1,1,1,1,1,1,
    2,2,2,2,1,1,1,1,1,
    2,2,1,1,2,2,1,1,1,
    4,0,2,0,2,0,2,1,0,
    8,0,4,0,2,0,2,1,0,
    8,0,2,0,4,0,2,1,0,
    16,0,4,0,4,0,2,1,0,
    4,4,2,2,2,2,1,1,1,
    16,0,4,0,4,0,4,1,0), byrow = T, ncol = 9)

  # If input is a pair of indivs, return the 9 coeffs as a numeric vector
  if(length(ids) == 2) {
    id1 = ids_int[1]; id2 = ids_int[2]
    RHS = c(
      1,
      2 * phi2(id1, id1, chromType = "x", mem = mem),
      2 * phi2(id2, id2, chromType = "x", mem = mem),
      4 * phi2(id1, id2, chromType = "x", mem = mem),
      8 * phi3(id1, id1, id2, chromType = "x", mem = mem),
      8 * phi3(id1, id2, id2, chromType = "x", mem = mem),
      16 * phi4(id1, id1, id2, id2, chromType = "x", mem = mem),
      4 * phi22(id1, id1, id2, id2, chromType = "x", mem = mem),
      16 * phi22(id1, id2, id1, id2, chromType = "x", mem = mem))

    j = solve(M9, RHS)

    # Set NA at undefined states (when males are involved)
    sex = getSex(x, c(id1, id2))
    if(sex[1] == 1 && sex[2] == 1)
      j[3:9] = NA
    if(sex[1] == 1 && sex[2] == 2)
      j[5:9] = NA
    if(sex[1] == 2 && sex[2] == 1)
      j[c(3:4,7:9)] = NA

    if(verbose)
      printCounts(mem)

    if(checkAnswer)
      compare_with_XIBD(x, ids, j)

    return(j)
  }
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


