#' Generalised kinship coefficients
#'
#' Computes single-locus generalised kinship coefficients of various kinds.
#' These are fundamental for computing identity coefficients (see
#' [identityCoefs()]), but are also interesting in their own right. Each
#' generalised kinship coefficient is defined as the probability of observing a
#' corresponding *generalised IBD pattern*, as defined and discussed in the
#' Details section below.
#'
#'
#' ## The starting point: standard kinship coefficients
#'
#' The classical kinship coefficient phi between two pedigree members A and B,
#' is the probability that two alleles sampled from A and B (one from each), at
#' a random autosomal locus, are identical by descent (IBD).
#'
#' In the language and notation to be introduced shortly, we would write `phi =
#' Pr[(A,B)]` where `(A,B)` is an *IBD pattern*.
#'
#' ## Generalised IBD patterns
#'
#' We define a *generalised IBD pattern* (GIP) to be a partition of a set of
#' alleles drawn from members of a pedigree, such that the alleles in each
#' subset are IBD. Each subset (also referred to as a *group* or a *block*) is
#' written as a collection of pedigree members (A, B, ...), with the
#' understanding that each member represents one of its alleles at the given
#' locus. A member may occur in multiple blocks, and also more than once within
#' a block.
#'
#' Additional requirements give rise to different flavours of GIPs (and their
#' corresponding coefficients):
#'
#' * `Distinct` (resp. `non-distinct`): alleles in different blocks are non-IBD
#' (resp. may be IBD)
#'
#' * `Deterministic` (resp. `random`): the parental origin (paternal or
#' maternal) of each allele is fixed (resp. unknown).
#'
#' We may say that a GIP is *partially* (rather than *fully*) deterministic if
#' the parental origin is fixed for some, but not all alleles involved.
#'
#'
#' ## Notational examples
#'
#' Our notation distinguishes the different types of patterns, as exemplified
#' below. Blocks are separated with "/" if they are distinct, and "&" otherwise.
#' Deterministically sampled alleles are suffixed by either ":p" (paternal) or
#' ":m" (maternal).
#'
#' * `(A, B) & (A, C)`: 4 alleles are sampled randomly; two from A, one from B
#' and one from C. The first from A is IBD with that from B, and the second from
#' A is IBD with that from C. All four alleles may be IBD. \[Random,
#' non-distinct\]
#'
#' * `(A, B) / (A, C)`: Same as the previous, but the two allele pairs must be
#' non-IBD. \[Random, distinct\]
#'
#' * `(A:p, C:p) / (C:m)`: The paternal alleles of A and C are IBD, and
#' different from the maternal allele of C. \[Deterministic, distinct\]
#'
#' * `(A, C:p) & (B, C:m)`: The paternal and maternal alleles of C are IBD with
#' random alleles of from A and B, respectively. The two pairs are not
#' necessarily different. \[Partially deterministic, non-distinct\]
#'
#' * `(A:p, A, A)`: Here we have just one group, specifying that paternal allele
#' of A is IBD with two other alleles sampled randomly from A. (If A is
#' non-inbred, this must have probability 1/4.) \[Partially deterministic,
#' single-block\]
#'
#' In the `kinpat()` constructor, deterministic sampling are indicated by naming
#' elements with "p" or "m". See Examples for how to create the above patterns.
#'
#'
#' ## Internal structure of `kinpat` objects
#'
#' (Note: This section is included only for completeness; `kinpat` objects
#' should not be directly manipulated by end users.)
#'
#' Internally, a GIP is stored as a list of integer vectors, each vector
#' giving the indices of pedigree members constituting an IBD block. In
#' addition, the object has three attributes:
#'
#' * `labs`: A character vector containing the names of all pedigree members
#'
#' * `deterministic`: A logical, which is TRUE if the pattern is (partially or
#' fully) deterministic
#'
#' * `distinct`: A logical.
#'
#' If `deterministic = TRUE`, a slight hack is used to preserve the additional
#' information. The parental origin is then encoded in the last digit (0 =
#' unknown; 1 = paternal; 2 = maternal), so that, e.g.:
#'
#' * 12 = the maternal origin of individual 1
#'
#' * 231 = the paternal allele of individual 23
#'
#' * 30 = a random allele of individual 3
#'
#'
#' ## A brief history of generalised kinship coefficients
#'
#' The notion of generalised kinship coefficients originated with Karigl (1981)
#' who used a selection of random, non-distinct patterns (in our terminology) to
#' compute identity coefficients.
#'
#' Weeks & Lange (1988), building on Karigl's work, defined random, distinct
#' patterns in full generality and gave an algorithm for computing the
#' corresponding coefficients.
#'
#' In a follow-up paper, Lange & Sinsheimer (1992) introduced partially
#' deterministic (distinct) patterns, and used these to compute detailed
#' identity coefficients.
#'
#' In another follow-up, Weeks et al. (1995) extended the work on random,
#' distinct patterns by Weeks & Lange (1988) to X-chromosomal loci.
#'
#' Garcia-Cortes (2015) gave an alternative algorithm for the detailed identity
#' coefficients, based on (in our terminology) fully deterministic, non-distinct
#' patterns.
#'
#'
#' ## Implemented algorithms
#'
#' The following are valid options for the `methods` parameters, and what they
#' implement.
#'
#' * `auto`: Chooses method automatically, based on the pattern type.
#'
#' * `K`: Karigl's algorithm for random, non-distinct patterns. Only a few cases
#' are supported, namely single groups up to length 4, and two groups of length
#' 2. (These were the ones considered by Karigl.) Extended to support
#' X-chromosomal patterns and inbred founders.
#'
#' * `WL`: Weeks & Lange's algorithm for random, distinct patterns of any size.
#' \[TODO: Include the extension to X by Weeks et al. (1995).\]
#'
#' * `LS`: Lange & Sinsheimer's algorithm for partially deterministic, distinct
#' patterns of any size. Does not support X, nor patterns involving inbred
#' founders.
#'
#' * `GC`: Garcia-Cortes' algorithm for fully deterministic, non-distinct
#' patterns. The current implementation only supports the patterns needed to
#' compute identity coefficients, namely single blocks and two blocks of length
#' 2. Extended to support X-chromosomal patterns and inbreed founders.
#'
#' @param x A `ped` object.
#' @param pattern A `kinpat` object, or a list of vectors to be passed onto
#'   [kinpat()]. Each vector should contain members of `x` constituting an IBD
#'   block. (See Details and Examples.)
#' @param distinct A logical indicating if different blocks are required to be
#'   non-IBD. Default: TRUE. (Irrelevant for single-block patterns.)
#' @param Xchrom A logical, by default FALSE.
#' @param method Either "auto", "K", "WL", "LS" or "GC".
#' @param verbose A logical, by default FALSE.
#' @param debug A logical, by default FALSE.
#' @param mem For internal use.
#' @param ... Further arguments.

#' @return `gKinship()` returns a single number, the probability of the given
#'   IBD pattern.
#'
#'   `kinpat()` returns an object of class `kinpat`. This is internally a list
#'   of integer vectors, with attributes `labs`, `deterministic` and `distinct`.
#'   (See also Details.)
#'
#' @seealso [kinship()], [identityCoefs()]
#'
#' @references
#'
#' * G. Karigl (1981). A recursive algorithm for the calculation of identity
#' coefficients. Ann. Hum. Genet.
#'
#' * D.E. Weeks & K. Lange (1988). The affected-pedigree-member method of
#' linkage analysis. Am. J. Hum. Genet
#'
#' * K. Lange & J.S. Sinsheimer (1992). Calculation of genetic identity
#' coefficients. Ann. Hum. Genet.
#'
#' * D.E. Weeks, T.I. Valappil, M. Schroeder, D.L. Brown (1995) An X-linked
#' version of the affected-pedigree-member method of linkage analysis. Hum
#' Hered.
#'
#' * L.A. García-Cortés (2015). A novel recursive algorithm for the calculation
#' of the detailed identity coefficients. Gen Sel Evol.
#'
#' @examples
#'
#' ### Trivial examples ###
#'
#' x = nuclearPed(father = "A", mother = "B", children = "C")
#'
#' # Random, distinct
#' kp1 = kinpat(x, list(c("A", "B"), c("A", "C")))
#' kp1
#'
#' # Random, non-distinct
#' kp2 = kinpat(x, list(c("A", "B"), c("A", "C")), distinct = FALSE)
#' kp2
#'
#' # Fully deterministic, distinct
#' kp3 = kinpat(x, list(c(p="A", p="C"), (m="C")))
#' kp3
#'
#' # Partially deterministic, non-distinct`
#' kp4 = kinpat(x, list(c("A", p="C"), c("B", m="C")), distinct = FALSE)
#' kp4
#'
#' # Partially deterministic, single block
#' kp5 = kinpat(x, list(c(p="A", "A", "A")))
#' kp5
#'
#' stopifnot(
#'   gKinship(x, kp1) == 0,      # (since A and B are unrelated)
#'   gKinship(x, kp2) == 0,      # (same as previous)
#'   gKinship(x, kp3) == 0.5,    # (only uncertainty is which allele A gave to C)
#'   gKinship(x, kp4) == 0.25,   # (distinct irrelevant)
#'   gKinship(x, kp5) == 0.25,   # (both random must hit the paternal)
#' )
#'
#'
#'
#' ### Kappa coefficients via generalised kinship ###
#'
#' # NB: Much less efficient than `kappaIBD()`; only for validation
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
#'
#' ### Detailed outputs and debugging ###
#'
#' x = fullSibMating(1)
#'
#' # Probability of sampling IBD alleles from 1, 5 and 6
#' p1 = kinpat(x, list(c(1,5,6)))
#' p1
#' gKinship(x, p1, method = "K", verbose = TRUE, debug = TRUE)
#' gKinship(x, p1, method = "WL", verbose = TRUE, debug = TRUE)
#'
#' # Probability that paternal of 5 is IBD with maternal of 6
#' p2 = kinpat(x, list(c(p=5, m=6)))
#' p2
#' gKinship(x, p2, method = "LS", verbose = TRUE, debug = TRUE)
#' gKinship(x, p2, method = "GC", verbose = TRUE, debug = TRUE)
#'
#' # Probability that paternal of 5 is *not* IBD with maternal of 6
#' p3 = kinpat(x, list(c(p=5), c(m=6)), distinct = TRUE)
#' p3
#' gKinship(x, p3, method = "LS", verbose = TRUE, debug = TRUE)
#'
#'
#' @export
gKinship = function(x, pattern, distinct = TRUE, Xchrom = FALSE,
                    method = c("auto", "K", "WL", "LS", "GC"),
                    verbose = FALSE, debug = FALSE, mem = NULL, ...) {

  method = match.arg(method)
  if(method == "auto")
    method = chooseIdentityMethod(x, pattern = pattern, Xchrom = Xchrom,
                                  detailed = isDeterministic(pattern))

  if(inherits(pattern, "kinpat"))
    pattern = kinpat2list(pattern)

  # Possibly add founder parents if method is LS
  addpar = NULL
  if(method == "LS") {
    ids = unlist(pattern)
    if(!is.null(names(ids)))
      addpar = ids[names(ids) != ""]
  }

  x = prepPed(x, addpar = addpar, Xchrom = Xchrom)

  # Now (re)make pattern object (point being: x might have changed)
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

