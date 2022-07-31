#' Expected LR of a pairwise kinship test
#'
#' Calculates the exact likelihood ratio of a pairwise kinship test,
#' implementing formulas of Egeland & Slooten (2016).
#'
#' @param x An hypothesised pedigree connecting two individuals.
#' @param true The true relationship between the two individuals.
#' @param ids A vector containing the names of the two individuals. Note: These
#'   must occur in both `x` and `true`.
#' @param L1 The number of alleles at the first locus.
#' @param L2 The number of alleles at the second locus, or NULL (default).
#' @param rho (If `L2` is not NULL.) A numeric vector of recombination
#'   fractions. Values outside the interval \eqn{[0, 0.5]} will raise an error.
#'
#' @return A single number, the expected LR.
#'
#' @references  Egeland, T. and Slooten, K. (2016). _The likelihood ratio as a
#'   random variable for linked markers in kinship analysis_. Int J Legal Med.
#'
#' @examples
#' #############################
#' # Fig. 2 of Egeland & Slooten
#' #############################
#'
#' rhos = seq(0, 0.5, length = 11)
#'
#' dat = cbind(
#'   Grand = ELR(linearPed(2), ids = c(1,5), L1 = 10, L2 = 30, rho = rhos),
#'   Half = ELR(halfSibPed(), ids = c(4,5), L1 = 10, L2 = 30, rho = rhos),
#'   Uncle = ELR(avuncularPed(), ids = c(3,6), L1 = 10, L2 = 30, rho = rhos))
#'
#' matplot(rhos, dat, type = "l",  lwd = 2, ylab = "E[LR]", ylim = c(0, 8))
#' legend("bottomleft", legend = colnames(dat), lty = 1:3, col = 1:3, lwd = 2)
#'
#'
#'
#' @export
ELR = function(x, true = x, ids = leaves(x), L1, L2 = NULL, rho = NULL) {

  if(!checkInt(L1))
    stop2("`L1` must be a positive integer: ", L1)

  if(linked <- !is.null(L2)) {
    if(!checkInt(L2))
      stop2("`L2` must be a positive integer: ", L2)
    if(!is.numeric(rho) || !all(rho >= 0 & rho <= 0.5))
      stop2("Recombination fraction(s) `rho` must be numbers in the interval [0, 0.5]")
  }

  # Matrix of Egeland & Slooten
  M = function(n) {
    matrix(c(1, 1, 1,
             1, (n+3)/4, (n+1)/2,
             1, (n+1)/2, n*(n+1)/2),
           ncol = 3)
  }

  ###  One locus
  if(!linked) {
    kap = kappaIBD(x, ids)
    kapTRUE = kappaIBD(true, ids)

    res = kap %*% M(L1) %*% t.default(kapTRUE)
    return(res)
  }

  ### Two linked loci: Formula (2.18) of Egeland & Slooten, 2016

  M1 = M(L1)
  M2 = M(L2)

  # Grid of 81 index combinations
  gr = expand.grid(rep(list(1:3), 4))

  res = sapply(rho, function(r) {
    K = twoLocusIBD(x, ids, rho = r)
    Ktrue = if(identical(x, true)) K else twoLocusIBD(true, ids, rho = r)
    terms = apply(gr, 1,
      function(a) K[a[1], a[2]] * Ktrue[a[3], a[4]] * M1[a[1], a[3]] * M2[a[2], a[4]])
    sum(terms)
  })

  res
}


checkInt = function(a) {
  length(a) == 1 && is.numeric(a) && a > 0 && as.integer(a) == a
}
