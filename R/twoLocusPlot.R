#' Two-locus coefficient plot
#'
#' Plot two-locus kinship or IBD coefficients as function of the recombination
#' rate.
#'
#' Each entry of `peds` must be a list with the following (named) entries:
#'
#' * ped: A ped object
#'
#' * ids: A pair of labels identifying two members of `ped`
#'
#' The `coeff` parameter must be either a character naming the coefficient to
#' compute, or a function. If a character, it must be one of the following
#' names: "inb", "kinship", "phi", "phi11", "k00", "k01", "k02", "k10", "k11", "k12",
#' "k20", "k21" or "k22".
#'
#' If `coeff` is a function, it must take three arguments named `ped`, `ids` and
#' `rho`, and produce a single number for each set of input data. See Examples.
#'
#' The first three are synonymous and indicate the two-locus kinship
#' coefficient. The remaining choices are two-locus IBD coefficients. (See
#' [twoLocusIBD()].)
#'
#' @param peds A list of lists. See details.
#' @param coeff A string identifying which coefficient to compute. See Details
#'   for legal values.
#' @param rseq A numeric vector of recombination rates. By default `seq(from =
#'   0, by = 0.5, length = 11)`.
#' @param xlab,ylab Axis labels.
#' @param col,lty,lwd Plotting parameters.
#' @param ... Further parameters passed on to [matplot()].
#'
#' @return NULL
#'
#' @examples
#'
#' ###############################
#' # Classic example of three relationships with equal one-locus coeffs
#' peds = list(
#'     GrandParent = list(ped = linearPed(2),    ids = c(1, 5)),
#'     HalfSib     = list(ped = halfSibPed(),    ids = c(4, 5)),
#'     Uncle       = list(ped = cousinPed(0, 1), ids = c(3, 6)))
#'
#' twoLocusPlot(peds, coeff = "kinship")
#' twoLocusPlot(peds, coeff = "k11")
#'
#' ###############################
#'
#' peds = list(
#'     PO = list(ped = nuclearPed(1), ids = c(1,3)),
#'     S  = list(ped = nuclearPed(2), ids = c(3,4)))
#'
#' twoLocusPlot(peds, coeff = "kinship")
#' twoLocusPlot(peds, coeff = "k11")
#'
#' ###############################
#'
#' ped1 = addChildren(halfSibPed(sex2 = 2), 4, 5, nch = 2)
#' ped2 = addChildren(addDaughter(nuclearPed(1), 3), 1, 5, nch = 2)
#' ped3 = addChildren(addDaughter(nuclearPed(2), 4), 3, 6, nch = 2)
#'
#' peds = list(
#'    `H-sibs` = list(ped = ped1, ids = leaves(ped1)),
#'    `G-sibs` = list(ped = ped2, ids = leaves(ped2)),
#'    `U-sibs` = list(ped = ped3, ids = leaves(ped3))
#' )
#' # plotPedList(peds)
#' twoLocusPlot(peds, coeff = "kinship")
#'
#' ################################
#'
#' ### Reproducing Fig 2 of Bishop & Williamson (1990)
#' ### This example illustrates `coeff` as a function.
#'
#' # The coefficient d11(rho) is the conditional probability of IBD = 1
#' # in the first locus, given IBD = 1 in the second.
#'
#' G = linearPed(2)
#' H = halfSibPed()
#' U = cousinPed(0, removal = 1)
#' FC = cousinPed(1)
#' FC1R = cousinPed(1, removal = 1)
#' SC = cousinPed(2)
#'
#' peds = list(
#'     GrandParent = list(ped = G,    ids = c(1, 5)),
#'     HalfSib     = list(ped = H,    ids = leaves(H)),
#'     Uncle       = list(ped = U,    ids = leaves(U)),
#'     FirstCous   = list(ped = FC,   ids = leaves(FC)),
#'     FirstCous1R = list(ped = FC1R, ids = leaves(FC1R)),
#'     SecondCous  = list(ped = SC,   ids = leaves(SC)))
#'
#'
#' d11 = function(ped, ids, rho) {
#'   twoLocusIBD(ped, ids, rho, coefs = "k11")/kappaIBD(ped, ids)[2]
#' }
#'
#' twoLocusPlot(peds, coeff = d11)
#'
#' @importFrom graphics legend matplot
#' @export
twoLocusPlot = function(peds, coeff = "k11", rseq = seq(0, 0.5, length = 11),
                        xlab = "Recombination rate",
                        ylab = NA, col = seq_along(peds), lty = 1, lwd = 1, ...) {

  if(!is.list(peds))
    stop2("Argument `peds` must be a list")

  if(is.ped(peds[[1]]))
    peds = list(peds)

  for(p in peds) {
    if(!is.list(p) || !setequal(names(p), c("ped", "ids")))
      stop2("Each entry of `peds` must be a list with names 'ped' and 'ids'")
    if(!is.ped(p$ped) || !all(p$ids %in% labels(p$ped)))
      stop2("Something is wrong with the input")
  }

  if(length(coeff) != 1)
    stop2("Argument `coeff` must have length 1: ", coeff)

  # Choose function
  if(is.function(coeff)) {
    myFUN = coeff
    if(is.na(ylab)) ylab = as.character(substitute(coeff))
  }
  else
    switch(coeff,
      inb = {
        myFUN = function(ped, ids, rho) twoLocusInbreeding(ped, id = ids, rho = rho)
        if(is.na(ylab)) ylab = "Two locus inbreeding"
      },
      kinship =, phi =, phi11 = {
        myFUN = function(ped, ids, rho) twoLocusKinship(ped, ids, rho = rho)
        if(is.na(ylab)) ylab = "Two locus kinship"
      },
      k00 =, k01 =, k02 =, k10 =, k11 =, k12 =, k20 =, k21 =, k22 = {
        myFUN = function(ped, ids, rho) twoLocusIBD(ped, ids, rho = rho, coefs = coeff)
        if(is.na(ylab)) ylab = coeff
      },
      stop2("Illegal value of `coeff`: ", coeff)
    )

  # Compute coefficients
  kvals = sapply(peds, function(x)
    sapply(rseq, function(r) myFUN(ped = x$ped, ids = x$ids, rho = r)))

  # Plot
  matplot(rseq, kvals, type = "l", xlab = xlab,
          ylab = ylab, col = col, lty = lty, lwd = lwd, ...)

  if(!is.null(names(peds)))
     legend("topright", names(peds), col = col, lty = lty, lwd = lwd)
}
