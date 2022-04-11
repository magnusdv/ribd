#' Generalised kinship coefficients
#'
#' Compute generalised single-locus kinship coefficients, as defined by Weeks &
#' Lange (1988), and also the further generalisation by Lange & Sinsheimer
#' (1992) allowing to specify if the alleles involved are paternally or
#' maternally inherited.
#'
#' @param x A `ped` object.
#' @param pattern A `kinpat` object, or a list to be passed onto [kinpat()].
#' @param Xchrom A logical, by default FALSE.
#' @param method Either "auto", "K", "WL", "LS" or "GC".
#' @param verbose A logical, by default FALSE.
#' @param mem For internal use.
#' @param debug A logical, by default FALSE.
#' @param ... Further arguments.
#'
#' @return A single probability.
#'
#' @examples
#'
#' x = nuclearPed(3)
#' kp = kinpat(x, list(c(1,1,1)))
#' gKinship(x, kp, method = "WL", debug = TRUE)
#'
#' # Deterministic
#' kp = kinpat(x, list(c(p=1,p=1), c(m=3, m=2)))
#' gKinship(x, kp, method = "GC", debug = TRUE)
#'
#'
#' ##### Kappa coefficients via generalised kinship ###
#'
#' # NB: Much less efficient than `kappaIBD()`; serves only as validation
#'
#' kappa_from_gk = function(x, ids, method = "WL") {
#'   fa1 = father(x, ids[1])
#'   fa2 = father(x, ids[2])
#'   mo1 = mother(x, ids[1])
#'   mo2 = mother(x, ids[2])
#'
#'   GK = function(...) gKinship(x, list(...), method = method)
#'
#'   k0 = GK(fa1, fa2, mo1, mo2)
#'   k1 = GK(c(fa1, fa2), mo1, mo2) + GK(c(fa1, mo2), fa2, mo1) +
#'        GK(c(mo1, fa2), fa1, mo2) + GK(c(mo1, mo2), fa1, fa2)
#'   k2 = GK(c(fa1, fa2), c(mo1, mo2)) + GK(c(fa1, mo2), c(mo1, fa2))
#'   c(k0, k1, k2)
#' }
#'
#' y1 = nuclearPed(2); ids = 3:4
#' stopifnot(kappa_from_gk(y1, ids) == kappaIBD(y1, ids))
#'
#' y2 = quadHalfFirstCousins(); ids = 9:10
#' stopifnot(kappa_from_gk(y2, ids) == kappaIBD(y2, ids))
#'
#' @export
gKinship = function(x, pattern, Xchrom = FALSE, mem = NULL, method = c("auto", "K", "WL", "LS", "GC"),
                    verbose = FALSE, debug = FALSE, ...) {

  method = match.arg(method)
  if(method == "auto")
    method = chooseIdentityMethod(x, unlist(pattern, use.names = FALSE), isDetailed(pattern), Xchrom)

  # Add founder parents if method is LS
  ids = if(method == "LS") unlist(pattern, use.names = FALSE) else NULL
  x = prepPed(x, addpar = ids, Xchrom = Xchrom)


  # If already kinpat object, reverse and remake (since x might have changed)
  if(inherits(pattern, "kinpat"))
    pattern = kinpat2list(pattern)
  kp = kinpat(x, pattern)

  if(is.null(mem))
    mem = memoIdentity(x, Xchrom = Xchrom, method = method, verbose = verbose, ...)

  res = switch(method,
               K = gKinship_Karigl(x, kp, Xchrom = Xchrom, mem = mem, debug = debug),
               WL = gKinship_WL(x, kp, Xchrom = Xchrom, mem = mem, debug = debug),
               LS = gKinship_LS(x, kp, Xchrom = Xchrom,mem = mem, debug = debug),
               GC = gKinship_GC(x, kp, Xchrom = Xchrom, mem = mem, debug = debug)
  )

  if(verbose)
    printMemInfo(mem)

  res
}

