#' IBD (kappa) coefficients
#'
#' Computes the three IBD coefficients summarising the relationship between two
#' non-inbred individuals. Both autosomal and X chromosomal versions are
#' implemented. The pedigree founders (other than the individuals in question)
#' are allowed to be inbred; see [pedtools::founderInbreeding()] for how to set
#' this up, and see Examples below.
#'
#' For non-inbred individuals a and b, their autosomal IBD coefficients
#' \eqn{(\kappa_0, \kappa_1, \kappa_2)} are defined as follows: \deqn{\kappa_i =
#' P(\text{a and b share exactly i alleles IBD at a random autosomal locus})}
#'
#' The autosomal kappa coefficients are computed from the kinship coefficients.
#' When a and b are both nonfounders, the following formulas hold:
#'
#' * \eqn{\kappa_2 = \varphi_{MM} \cdot \varphi_{FF} + \varphi_{MF} \cdot\varphi_{FM}}
#'
#' * \eqn{\kappa_1 = 4 \varphi_{ab} - 2 \kappa_2}
#'
#' * \eqn{\kappa_0 = 1 - \kappa_1 - \kappa_2}
#'
#' Here \eqn{\varphi_{MF}} denotes the kinship coefficient between the
#' **m**other of a and the **f**ather of b, etc. If either a or b is a founder,
#' then \eqn{\kappa_2 = 0}, while the other two formulas remain as before.
#'
#' The X-chromosomal IBD coefficients are defined similarly to the autosomal
#' case. Here \eqn{\kappa_2} is undefined when one or both individuals are male,
#' which greatly simplifies the calculations when males are involved. The
#' formulas are (with \eqn{\varphi_{ab}} now referring to the X-chromosomal
#' kinship coefficient):
#'
#' * Both male: \eqn{(\kappa_0, \kappa_1, \kappa_2) = (1-\varphi_{ab}, \varphi_{ab}, \text{NA})}
#'
#' * One male, one female: \eqn{(\kappa_0, \kappa_1, \kappa_2) = (1-2 \varphi_{ab},
#' 2 \varphi_{ab}, \text{NA})}
#'
#' * Two females: Similar formulas as in the autosomal case.
#'
#' @param x A pedigree in the form of a `ped` object (or a list of such).
#' @param ids A character (or coercible to character) containing ID labels of
#'   two or more pedigree members.
#' @param inbredAction An integer telling the program what to do if either of
#'   the `ids` individuals are inbred. Possible values are: 0 = do nothing; 1 =
#'   print a warning message (default); 2 = raise an error. In the first two
#'   cases the coefficients are reported as `NA`.
#' @param simplify Simplify the output (to a numeric of length 3) if `ids` has
#'   length 2. Default: TRUE.
#' @param acrossComps A logical indicating if pairs of individuals in different
#'   components should be included. Default: TRUE.
#' @param Xchrom A logical, indicating if the autosomal (default) or
#'   X-chromosomal kappa coefficients should be computed.
#'
#' @return If `ids` has length 2 and `simplify = TRUE`: A numeric vector of
#'   length 3: \eqn{(\kappa_0, \kappa_1, \kappa_2)}.
#'
#'   Otherwise: A data frame with one row for each pair of individuals, and 5
#'   columns. The first two columns contain the ID labels, and columns 3-5
#'   contain the IBD coefficients.
#'
#'   Kappa coefficients of inbred individuals (meaning X-inbred females if
#'   `Xchrom = T`) are reported as NA, unless `inbredAction = 2`.
#'
#'   The X-chromosomal \eqn{\kappa_2} is NA whenever at least one of the two
#'   individuals is male.
#'
#' @seealso [kinship()], [identityCoefs()] for coefficients allowing inbreeding,
#'   [showInTriangle()] for plotting kappa coefficients in the IBD triangle.
#'
#' @examples
#' ### Siblings
#' x = nuclearPed(2)
#' kappaIBD(x)
#'
#' k = kappaIBD(x, 3:4)
#' stopifnot(identical(k, c(.25, .5, .25)))
#'
#' ### Quad half first cousins
#' x = quadHalfFirstCousins()
#' k = kappaIBD(x, ids = leaves(x))
#' stopifnot(identical(k, c(17/32, 14/32, 1/32)))
#'
#' ### Paternal half brothers with 100% inbred father
#' # Genetically indistinguishable from an (outbred) father-son relationship
#' x = halfSibPed() |> setFounderInbreeding(ids = 2, value = 1)
#' plot(x, hatched = 4:5)
#'
#' k = kappaIBD(x, 4:5)
#' stopifnot(identical(k, c(0, 1, 0)))
#'
#' ### X-chromosomal kappa
#' y = nuclearPed(2, sex = 2)
#' kappaIBD(y, Xchrom = TRUE)
#'
#' @export
kappaIBD = function(x, ids = labels(x), inbredAction = 1, simplify = TRUE,
                    acrossComps = TRUE, Xchrom = FALSE) {

  if(is.pedList(x)) {
    compNr = getComponent(x, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
    compNr = unique.default(compNr)

    if(length(compNr) == 1)
      return(kappaIBD(x[[compNr]], ids, inbredAction = inbredAction, simplify = simplify, Xchrom = Xchrom))

    x = x[compNr]
    nPed = length(x)
    idsComp = lapply(x, function(comp) intersect(comp$ID, ids))

    # Within-component coefficients
    kapComp = lapply(which(lengths(idsComp) > 1), function(i)
      kappaIBD(x[[i]], idsComp[[i]], inbredAction = inbredAction, simplify = FALSE, Xchrom = Xchrom))
    res = do.call(rbind, kapComp)

    # Between-components
    if(acrossComps) {
      for(i in seq_len(nPed - 1)) for(j in seq(i+1, nPed)) for(a in idsComp[[i]])
        res = rbind(res, data.frame(id1 = a, id2 = idsComp[[j]], kappa0 = 1, kappa1 = 0, kappa2 = 0, stringsAsFactors = FALSE))
    }

    # Simplify output for a single pair
    if(simplify && length(ids) == 2)
      res = as.numeric(res[1, 3:5])

    return(res)
  }

  # If X, delegate to X version
  if(Xchrom)
    return(.kappaX(x, ids, inbredAction = inbredAction, simplify = simplify))

  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  labs = labels(x)
  ids = as.character(ids)
  if(!all(ids %in% labs))
    stop2("Unknown ID label: ", setdiff(ids, labs))
  if(length(ids) < 2)
    stop2("At least two ID labels must be indicated")
  if(dup <- anyDuplicated.default(ids))
    stop2("Duplicated ID label: ", ids[dup])

  KIN = kinship(x)
  INB = 2*diag(KIN) - 1 # inbreeding coeffs

  isInbred = INB[ids] > .Machine$double.eps
  if(any(isInbred) && inbredAction > 0) {
    msg = paste0(c(sprintf(" Individual '%s' is inbred (f = %g)", ids[isInbred], INB[ids[isInbred]]),
                  "Kappa coefficients are only defined for non-inbred individuals."), collapse = "\n")
    switch(inbredAction, message(msg), stop2(msg))
  }

  # Build result data frame
  pairs = .comb2(ids, vec = TRUE)
  res = data.frame(id1 = pairs[, 1], id2 = pairs[, 2],
                   kappa0 = NA_real_, kappa1 = NA_real_, kappa2 = NA_real_,
                   stringsAsFactors = FALSE)

  # Rows that needs computing
  eps = .Machine$double.eps
  noninbred_rows = INB[res$id1] < eps & INB[res$id2] < eps

  fous = founders(x)
  fou_rows = res$id1 %in% fous | res$id2 %in% fous

  # Noninbred pairs involving founder(s)
  fou_noninb_rows = fou_rows & noninbred_rows
  if(any(fou_noninb_rows)) {
    phi_fou = KIN[pairs[fou_noninb_rows, , drop = FALSE]]
    res[fou_noninb_rows, 3:5] = cbind(1 - 4*phi_fou, 4*phi_fou, 0)
  }

  # Noninbred nonfounders
  if(any(nn_rows <- !fou_rows & noninbred_rows)) {
    id1.nn = res$id1[nn_rows]
    id2.nn = res$id2[nn_rows]
    F1 = father(x, id1.nn)
    M1 = mother(x, id1.nn)
    F2 = father(x, id2.nn)
    M2 = mother(x, id2.nn)

    k2 = KIN[cbind(F1, F2)]*KIN[cbind(M1, M2)] + KIN[cbind(F1, M2)]*KIN[cbind(M1, F2)]
    k1 = 4*KIN[pairs[nn_rows, , drop = FALSE]] - 2*k2
    k0 = 1 - k1 - k2
    res[nn_rows, 3:5] = cbind(k0, k1, k2)
  }

  # Simplify output for a single pair
  if(simplify && length(ids) == 2)
    res = as.numeric(res[1, 3:5])

  res
}


.kappaX = function(x, ids, inbredAction = 1, simplify = TRUE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  labs = labels(x)
  ids = as.character(ids)
  if(!all(ids %in% labs))
    stop2("Unknown ID label: ", setdiff(ids, labs))
  if(length(ids) < 2)
    stop2("At least two ID labels must be indicated")
  if(dup <- anyDuplicated.default(ids))
    stop2("Duplicated ID label: ", ids[dup])

  SEX = getSex(x, named = TRUE)
  KIN = kinship(x, Xchrom = TRUE)
  INB = 2*diag(KIN) - 1 # inbreeding coeffs

  isInbred = SEX[ids] == 2 & INB[ids] > .Machine$double.eps
  if(any(isInbred) && inbredAction > 0) {
    msg = paste0(c(sprintf("* Individual '%s' is X-inbred (f = %g)", ids[isInbred], INB[ids[isInbred]]),
                   "X-chromosomal kappa coefficients are only defined for pairs not involving X-inbred females."), collapse = "\n")
    switch(inbredAction, message(msg), stop2(msg))
  }

  # Build result data frame
  pairs = .comb2(ids, vec = TRUE)
  res = data.frame(id1 = pairs[, 1], id2 = pairs[, 2],
                   kappa0 = NA_real_, kappa1 = NA_real_, kappa2 = NA_real_,
                   stringsAsFactors = FALSE)

  # Inbred rows should remain NA
  inbred1 = SEX[res$id1] == 2 & INB[res$id1] > .Machine$double.eps
  inbred2 = SEX[res$id2] == 2 & INB[res$id2] > .Machine$double.eps
  inbred_rows = inbred1 | inbred2


  ### Male-male
  mm_rows = SEX[res$id1] == 1 & SEX[res$id2] == 1
  if(any(mm_rows)) {
    phi_mm = KIN[pairs[mm_rows, , drop = FALSE]]
    res[mm_rows, 3:5] = cbind(1 - phi_mm, phi_mm, NA_real_)
  }

  ### Male-female (or vice versa)
  mf_rows = (SEX[res$id1] == 1 & SEX[res$id2] == 2) | (SEX[res$id1] == 2 & SEX[res$id2] == 1)
  mf_rows = mf_rows & !inbred_rows

  if(any(mf_rows)) {
    phi_mf = KIN[pairs[mf_rows, , drop = FALSE]]
    res[mf_rows, 3:5] = cbind(1 - 2*phi_mf, 2*phi_mf, NA_real_)
  }

  #### Female-female
  ff_rows = (SEX[res$id1] == 2 & SEX[res$id2] == 2) & !inbred_rows

  fous = founders(x)
  founder_rows = res$id1 %in% fous | res$id2 %in% fous

  # Pairs involving at least one founder
  ff_fou_rows = ff_rows & founder_rows
  if(any(ff_fou_rows)) {
    phi_ff_fou = KIN[pairs[ff_fou_rows, , drop = FALSE]]
    res[ff_fou_rows, 3:5] = cbind(1 - 4*phi_ff_fou, 4*phi_ff_fou, 0)
  }

  # Both nonfounders
  ff_nonfou_rows = ff_rows & !founder_rows
  if(any(ff_nonfou_rows)) {
    id1_ff_nonfou = res$id1[ff_nonfou_rows]
    id2_ff_nonfou = res$id2[ff_nonfou_rows]
    F1 = father(x, id1_ff_nonfou)
    M1 = mother(x, id1_ff_nonfou)
    F2 = father(x, id2_ff_nonfou)
    M2 = mother(x, id2_ff_nonfou)

    k2 = KIN[cbind(F1, F2)]*KIN[cbind(M1, M2)] + KIN[cbind(F1, M2)]*KIN[cbind(M1, F2)]
    k1 = 4*KIN[pairs[ff_nonfou_rows, , drop = FALSE]] - 2*k2
    k0 = 1 - k1 - k2
    res[ff_nonfou_rows, 3:5] = cbind(k0, k1, k2)
  }

  # Simplify output for a single pair
  if(simplify && length(ids) == 2)
    res = as.numeric(res[1, 3:5])

  res
}

