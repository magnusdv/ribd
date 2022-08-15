#' Variance of realised relatedness coefficients
#'
#' Compute the variance of realised relatedness coefficients, by doubly
#' integrating the corresponding two-locus coefficients.
#'
#' The double integral method was used by Guo to compute the variation in genome
#' sharing (Guo 1995, see also Thompson 2013). The method extends directly to
#' other coefficients. The implementation here supports Cotterman's kappa
#' coefficients (of noninbred individuals), and Jacquard's condensed identity
#' coefficients.
#'
#' The current implementation is based on stats::integrate and can probably be
#' optimised in several ways.
#'
#' @param x A ped object.
#' @param ids A vector naming two members of `x`.
#' @param coef A single character: one of "k0", "k1", "k2" in `kappaVariance`,
#'   and one of "D1", "D2", ... "D9" in `deltaVariance`.
#' @param M A positive number; the chromosome length in Morgan.
#'
#' @return A positive number.
#'
#' @references
#'
#' * Guo (1995) _Proportion of genome shared identical by descent by relatives:
#' concept, computation, and applications_. Am J Hum Genet.
#'
#' * Hill & Weir (2011). _Variation in actual relationship as a consequence of
#' Mendelian sampling and linkage_. Genet Res.
#'
#' * Thompson (2013). _Identity by Descent: Variation in Meiosis, Across
#' Genomes, and in Populations_. Genetics.
#'
#' @examples
#' ###################################
#' ### Box 1 of Hill & Weir (2011) ###
#' ###################################
#'
#' # Eq. 4b of Hill & Weir
#' phi = function(n, l) {
#'  1/(2*l^2) * (1/4)^n * sum(sapply(1:n, function(r)
#'    choose(n, r) * (2*r*l - 1 + exp(-2*r*l))/r^2))
#' }
#'
#' # Chromosome of 1 Morgan
#' M = 1
#'
#' ### Full sibs ###
#'
#' x = nuclearPed(2)
#' kappaVariance(x, ids = 3:4, coef = "k2", M = M)
#'
#' # Hill & Weir (Box 1)
#' 16*phi(4,M) - 16*phi(3,M) + 8*phi(2,M) - 2*phi(1,M)
#'
#'
#' ### Double first cousins ###
#'
#' \dontrun{
#' dfc = doubleFirstCousins()
#' kappaVariance(dfc, coef = "k0", M = M)
#' kappaVariance(dfc, coef = "k1", M = M)
#' kappaVariance(dfc, coef = "k2", M = M)
#'
#' # Hill & Weir, Box 1
#' var_k2 = 64*phi(8,M) - 64*phi(7,M) + 40*phi(6,M) - 20*phi(5,M) +
#'  33/4*phi(4,M) - 5/2*phi(3,M) + 5/8*phi(2,M)-1/8*phi(1,M)
#' var_k1 = 4*var_k2
#' var_k0 = var_k2 + 2 * (4*phi(4,M) - 2*phi(3,M) + 3/4*phi(2,M) - 1/4*phi(1,M))
#'
#' var_k0
#' var_k1
#' var_k2
#'
#' }
#'
#' @importFrom stats integrate qchisq var
#' @export
kappaVariance = function(x, ids = leaves(x), coef, M = 1) {
  coef = match.arg(coef, c("k0", "k1", "k2"))
  coef2 = switch(coef, k0 = "k00", k1 = "k11", k2 = "k22")

  # Double integral, as in Thompson 2013
  E2 = 2/M^2 *
    integrate(Vectorize(function(a) {
      integrate(Vectorize(function(b) {
        rho = 0.5 * (1 - exp(-2*(a-b)))
        twoLocusIBD(x, ids, rho = rho, coefs = coef2)
      }), 0, a)$value
    }), 0, M)$value

  # Single locus coefficient
  kappa1L = kappaIBD(x, ids)[match(coef, c("k0", "k1", "k2"))]

  # Var = E(X^2) - E(X)^2
  var = E2 - kappa1L^2
  list(mean = kappa1L, var = var)
}


#' @export
#' @rdname kappaVariance
deltaVariance = function(x, ids = leaves(x), coef, M = 1) {
  coef = match.arg(coef, paste0("D", 1:9))

  # Double integral, as in Thompson 2013
  E2 = 2/M^2 *
    integrate(Vectorize(function(a) {
      integrate(Vectorize(function(b) {
        rho = 0.5 * (1 - exp(-2*(a-b)))
        twoLocusIdentity(x, ids, rho = rho)[[coef, coef]]
      }), 0, a)$value
    }), 0, M)$value

  # Single locus coefficient
  delta1L = identityCoefs(x, ids)[match(coef, paste0("D", 1:9))]

  # Var = E(X^2) - E(X)^2
  var = E2 - delta1L^2
  list(mean = delta1L, var = var)
}



# kappaVarSim = function(x, ids = leaves(x), coef = "k1", M = 1, nsim = 1000) {
#   map = uniformMap(M = M)
#   sims = ibdsim(x, N = nsim, ids = ids, map = map, model = "haldane", verbose = F)
#   kap = realisedKappa(sims)$perSimulation[, coef]
#   list(mean = mean(kap),
#        var = var(kap),
#        CIvar = (nsim - 1)*var(kap) / qchisq(c(.975,.025), nsim-1))
# }

