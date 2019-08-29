#' Two-locus identity coefficients
#'
#' Computes the 9*9 matrix of two-locus condensed identity coefficients of a
#' pair of pedigree members, for a given recombination rate.
#'
#' Let A, B be two pedigree members, and L1, L2 two loci with a given
#' recombination rate \eqn{\rho}. The two-locus identity coefficients
#' \eqn{\Delta_{i,j}(\rho)}{\Delta_ij(\rho)}, for \eqn{1 \le i,j \le 9} are
#' defined as the probability that the identity state of the alleles of A and B
#' are \eqn{\Sigma_i} at L1 and \eqn{\Sigma_j} at L2 simultaneously. (The
#' ordering of the 9 states follows Jacquard (1974).)
#'
#' The method of computation depends on whether the relationship between A and B
#' is rectilineal or not. Only non-rectilineal relationships are implemented for
#' now.
#'
#' @param x A pedigree in the form of a [`pedtools::ped`] object.
#' @param ids A character (or coercible to character) containing ID labels of
#'   two pedigree members.
#' @param rho A number in the interval \eqn{[0, 0.5]}; the recombination rate
#'   between the two loci.
#' @param coefs A character indicating which coefficient(s) to compute. A subset
#'   of `c('d00', 'd01', ..., 'd99')`. By default, all coefficients are
#'   computed.
#' @param detailed A logical, indicating whether the condensed (default) or
#'   detailed coefficients should be returned.
#' @param verbose A logical.
#'
#' @return By default, a symmetric 9*9 matrix containing the two-locus condensed
#'   identity coefficients \eqn{\Delta_{i,j}}{\Delta_ij}.
#'
#'   If either `coefs` is explicitly given (i.e., not NULL), or `detailed =
#'   TRUE`, the computed coefficients are returned as a named vector.
#'
#' @seealso [twoLocusIBD()]
#'
#' @examples
#' ### Full sibs ###
#' x = nuclearPed(2)
#' kapp = twoLocusIBD(x, ids = 3:4, rho = 0.25)
#' jacq = twoLocusIdentity(x, ids = 3:4, rho = 0.25)
#' stopifnot(all.equal(jacq[9:7,9:7], kapp, check.attributes = FALSE))
#'
#' ### Full sib mating ###
#' x = fullSibMating(1)
#' j = condensedIdentity(x, ids = 5:6)
#' j2 = twoLocusIdentity(x, ids = 5:6, rho = 0.25)
#' stopifnot(identical(unname(rowSums(j2)), j))
#'
#'
#' @export
twoLocusIdentity = function(x, ids, rho, coefs = NULL, detailed = F, verbose = F) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must have length exactly 2")

  if(!is.null(coefs)) stop2("Argument `coefs` is not implemented yet")
  if(detailed) stop2("Argument `detailed` is not implemented yet")

  # One-locus identity coefficients
  # j1 = condensedIdentity(x, ids)

  # Enforce parents to precede their children
  if(!hasParentsBeforeChildren(x))
    x = parentsBeforeChildren(x)

  x = foundersFirst(x)

  ids = internalID(x, ids)
  a = ids[1]
  b = ids[2]

  # Setup memoisation
  mem = initialiseTwoLocusMemo(x, rho, counters = c("i", "ilook", "ir", "b0", "b1", "b2", "b3"))

  # Output format
  outputMatrix = is.null(coefs) && !detailed

  # Coefficient selection
  #allcoefs = c("k00", "k01", "k02", "k10", "k11", "k12", "k20", "k21", "k22")
  #if(any(!coefs %in% allcoefs))
  #  stop2("Invalid coefficient specified: ", setdiff(coefs, allcoefs), "\nSee ?twoLocusIdentity for valid choices.")
  #if(is.null(coefs))
  #  coefs = allcoefs

  # Rectilineal
  rect_a = a %in% ancestors(x, b, internal = T)
  rect_b = b %in% ancestors(x, a, internal = T)
  if(rect_a || rect_b)
    stop2("Rectilineal relationships are not implemented yet")

  if(!rect_a && !rect_b) {
    RES = twoLocusIdentity_lateral(x, ids, rho, mem = mem, coefs = coefs, detailed = detailed, verbose = verbose)
  }

  ### Print summary
  if(verbose)
    printCounts2(mem)

  # Output
  if(outputMatrix) {
    dim(RES) = c(9, 9)
    dimnames(RES) = list(paste0("D", 1:9), paste0("D", 1:9))
  }

  RES
}


twoLocusIdentity_lateral = function(x, ids, rho, mem = NULL, coefs, detailed = F, verbose = F) {

  if(is.null(mem)) {
    # Enforce parents to precede their children
    if(!hasParentsBeforeChildren(x))
      x = parentsBeforeChildren(x)

    x = foundersFirst(x)

    # Setup memoisation
    mem = initialiseTwoLocusMemo(x, rho, counters = c("i", "ilook", "ir", "b1", "b2", "b3"))
  }

  idsi = internalID(x, ids)
  a = idsi[1]; b = idsi[2]

  # If unrelated, return 0 matrix
  if(mem$k1[a, b] == 0) {
    RES = matrix(0, ncol = 9, nrow = 9, dimnames = list(paste0("D", 1:9), paste0("D", 1:9)))
    if(mem$isCompletelyInbred[a] && mem$isCompletelyInbred[b])
      RES[2,2] = 1
    else if(mem$isCompletelyInbred[a])
      RES[4,4] = 1
    else if(mem$isCompletelyInbred[b])
      RES[6,6] = 1
    else
      RES[9,9] = 1

    return(RES)
  }

  # By now, none should be founders
  if(any(ids %in% founders(x)))
    stop2("The lateral method cannot be used for these individuals. Rectilineal relationship?")

  # Parents (internal indices)
  f = father(x, a, internal = T)
  g = father(x, b, internal = T)
  m = mother(x, a, internal = T)
  n = mother(x, b, internal = T)

  # Meiosis indicators
  # Parental meioses must be separated in case of selfing
  f.mei = sprintf("%s>%s", f, 100*a + 1)
  m.mei = sprintf("%s>%s", m, 100*a + 2)

  g.mei = sprintf("%s>%s", g, 100*b + 1)
  n.mei = sprintf("%s>%s", n, 100*b + 2)

    # Detailed single-locus identity states, described as generalised kinship patterns
  S = c(
    sprintf("%s = %s = %s = %s", f.mei, g.mei, m.mei, n.mei), # 1
    sprintf("%s = %s = %s , %s", f.mei, g.mei, m.mei, n.mei), # 2
    sprintf("%s = %s = %s , %s", f.mei, n.mei, m.mei, g.mei), # 3
    sprintf("%s = %s = %s , %s", f.mei, g.mei, n.mei, m.mei), # 4
    sprintf("%s = %s = %s , %s", m.mei, g.mei, n.mei, f.mei), # 5
    sprintf("%s = %s , %s = %s", f.mei, m.mei, g.mei, n.mei), # 6
    sprintf("%s = %s , %s , %s", f.mei, m.mei, g.mei, n.mei), # 7
    sprintf("%s = %s , %s , %s", g.mei, n.mei, f.mei, m.mei), # 8

    sprintf("%s = %s , %s = %s", f.mei, g.mei, m.mei, n.mei), # 9
    sprintf("%s = %s , %s , %s", f.mei, g.mei, m.mei, n.mei), # 10
    sprintf("%s = %s , %s , %s", m.mei, n.mei, f.mei, g.mei), # 11
    sprintf("%s = %s , %s = %s", f.mei, n.mei, m.mei, g.mei), # 12
    sprintf("%s = %s , %s , %s", f.mei, n.mei, m.mei, g.mei), # 13
    sprintf("%s = %s , %s , %s", m.mei, g.mei, f.mei, n.mei), # 14
    sprintf("%s , %s , %s , %s", f.mei, g.mei, m.mei, n.mei)  # 15
  )

  # Which detailed states belong to each condensed state (Jacquard ordering!)
  I = list(
    1,
    6,
    2:3,
    7,
    4:5,
    8,
    c(9,12),
    c(10,11,13,14),
    15
  )

  # Helper function for computing generalised two-locus coefs
  .Phi = function(loc1, loc2) {
    kin = kin2L(x, loc1, loc2, internal = T)
    genKin2L(kin, mem, indent = NA)
  }

  ### Here starts the actual work! ###

  res = matrix(0, nrow = 9, ncol = 9)
  for(i in 1:9) for(j in i:9) {
    for(u in I[[i]]) for(v in I[[j]])
      res[i,j] = res[i,j] + .Phi(S[u], S[v])

    res[j,i] = res[i,j]
  }

  res
}

