#' Variance of realised relatedness coefficients
#'
#' Compute the variance of realised relatedness coefficients, by doubly
#' integrating the corresponding two-locus coefficients.
#'
#' The double integral method was used by Guo to compute the variation in
#' proportion of the genome shared IBD (Guo 1995, see also Thompson 2013). The
#' method extends directly to other coefficients. The implementation here
#' supports Cotterman's kappa coefficients (of noninbred individuals), and
#' Jacquard's condensed identity coefficients.
#'
#' This function is a bare-bones implementation of the double integral method,
#' based on `stats::integrate`, and can probably be optimised in various ways.
#'
#' The `coeff` parameter must be either a character naming the coefficient to
#' compute, or a function. If a character, it must be one of the following
#' names:
#'
#' * "inb" (inbreeding coefficient)
#'
#' * "kinship", "phi" (synonyms for the kinship coefficient)
#'
#' * "k0", "k1", "k2" (kappa coefficients of noninbred individuals)
#'
#' * "D1", "D2", ... "D9" (condensed identity coefficients)
#'
#' @param x A ped object.
#' @param ids A vector naming two members of `x`.
#' @param coeff A string naming a coefficient for which the variance is to be
#'   computed. See Details for legal values.
#' @param L A positive number; the chromosome length in Morgan.
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
#' L = 1
#'
#' ### Full sibs ###
#'
#' x = nuclearPed(2)
#' realisedIbdVariance(x, ids = 3:4, coeff = "k2", L = L)
#'
#' # Hill & Weir (Box 1)
#' 16*phi(4,L) - 16*phi(3,L) + 8*phi(2,L) - 2*phi(1,L)
#'
#'
#' ### Double first cousins ###
#'
#' \dontrun{
#' dfc = doubleFirstCousins()
#'
#' # Runtime ~1 min
#' realisedIbdVariance(dfc, coeff = "k0", L = L)
#' realisedIbdVariance(dfc, coeff = "k1", L = L)
#' realisedIbdVariance(dfc, coeff = "k2", L = L)
#'
#' # Hill & Weir, Box 1
#' var_k2 = 64*phi(8,L) - 64*phi(7,L) + 40*phi(6,L) - 20*phi(5,L) +
#'  33/4*phi(4,L) - 5/2*phi(3,L) + 5/8*phi(2,L)-1/8*phi(1,L)
#' var_k1 = 4*var_k2
#' var_k0 = var_k2 + 2 * (4*phi(4,L) - 2*phi(3,L) + 3/4*phi(2,L) - 1/4*phi(1,L))
#'
#' var_k0
#' var_k1
#' var_k2
#'
#' }
#'
#' @importFrom stats integrate qchisq var
#' @export
realisedIbdVariance = function(x, ids = leaves(x), coeff, L = 1) {

  # Choose function
  switch(coeff,
         inb = {
           mu = inbreeding(x, ids)
           twoLocFun = function(rho) twoLocusInbreeding(x, ids, rho = rho)
         },
         kinship =, phi = {
           mu = kinship(x, ids)
           twoLocFun = function(rho) twoLocusKinship(x, ids, rho = rho)
         },
         R = { # "actual relationship" (Hill & Weir, 2011)
           kap = kappaIBD(x, ids) # raise error if inbreding
           mu = kap[2]/2 + kap[3]  # = 2*phi
           twoLocFun = function(rho) {
             K = twoLocusIBD(x, ids, rho = rho)
             K[2,2]/4 + K[2,3] + K[3,3]
           }
         },
         k0 = {
           mu = kappaIBD(x, ids)[1]
           twoLocFun = function(rho) twoLocusIBD(x, ids, rho = rho, coefs = "k00")
         },
         k1 = {
           mu = kappaIBD(x, ids)[2]
           twoLocFun = function(rho) twoLocusIBD(x, ids, rho = rho, coefs = "k11")
         },
         k2 = {
           mu = kappaIBD(x, ids)[3]
           twoLocFun = function(rho) twoLocusIBD(x, ids, rho = rho, coefs = "k22")
         },
         D1 =, D2 =, D3 =, D4 =, D5 =, D6 =, D7 =, D8 =, D9 = {
           mu = identityCoefs(x, ids)[match(coeff, paste0("D", 1:9))]
           twoLocFun = function(rho) twoLocusIdentity(x, ids, rho = rho)[[coeff, coeff]]
         },
         stop2("Illegal value of `coeff`: ", coeff)
  )

  # Double integral, as in Thompson 2013 (but integrating over triangle)
  E2 = 2/L^2 *
    integrate(Vectorize(function(a) {
      integrate(Vectorize(function(b) {
        rho = 0.5 * (1 - exp(-2*(a-b)))
        twoLocFun(rho)
      }), 0, a)$value
    }), 0, L)$value

  # Var(X) = E(X^2) - E(X)^2
  list(mean = mu, var = E2 - mu^2)
}



# kappaVarSim = function(x, ids = leaves(x), coef = "k1", L = 1, nsim = 1000) {
#   map = uniformMap(L = L)
#   sims = ibdsim(x, N = nsim, ids = ids, map = map, model = "haldane", verbose = F)
#   kap = realisedKappa(sims)$perSimulation[, coef]
#   list(mean = mean(kap),
#        var = var(kap),
#        CIvar = (nsim - 1)*var(kap) / qchisq(c(.975,.025), nsim-1))
# }

