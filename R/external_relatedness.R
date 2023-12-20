#' Relatedness coefficients by other programs
#'
#' Wrappers for functions in other packages or external programs.
#'
#' `kinship2_kinship()` and `kinship2_inbreeding()` both wrap
#' [kinship2::kinship()].
#'
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#' @param ids A integer vector of length 2.
#' @param Xchrom A logical, indicating if the autosomal (default) or
#'   X-chromosomal coefficients should be computed.
#'
#' @return
#'
#' For `kinship2_inbreeding()`, a numerical vector with inbreeding coefficients,
#' named with ID labels.
#'
#' For `kinship2_kinship()`, either a single numeric (if `ids` is a pair of
#' pedigree members) or the whole kinship matrix, with the ID labels as
#' dimnames.
#'
#' @seealso [kinship2::kinship()]
#'
#' @examples
#' # A random pedigree with 7 individuals
#' p = randomPed(n = 7, seed = 123)
#'
#' ### Kinship matrix
#'
#' # Autosomal: Check that ribd agrees with kinship2
#' stopifnot(identical(
#'   kinship(p),          # ribd
#'   kinship2_kinship(p)  # kinship2
#' ))
#'
#' # X chromosomal kinship
#' stopifnot(identical(
#'   kinship(p, Xchrom = TRUE),          # ribd
#'   kinship2_kinship(p, Xchrom = TRUE)  # kinship2
#' ))
#'
#'
#' ### Inbreeding coefficients
#'
#' # Autosomal
#' stopifnot(identical(
#'   inbreeding(p),          # ribd
#'   kinship2_inbreeding(p)  # kinship2
#' ))
#'
#' # X chromosomal
#' stopifnot(identical(
#'   inbreeding(p, Xchrom = TRUE),          # ribd
#'   kinship2_inbreeding(p, Xchrom = TRUE)  # kinship2
#' ))
#'
#'
#' @name external_coefs
NULL

#' @rdname external_coefs
#' @importFrom kinship2 kinship
#' @export
kinship2_kinship = function(x, ids = NULL, Xchrom = FALSE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  kin.matrix = kinship2::kinship(1:pedsize(x), dadid = x$FIDX, momid = x$MIDX,
                                 sex = x$SEX, chrtype = if(Xchrom) "x" else "autosome")
  labs = labels(x)
  dimnames(kin.matrix) = list(labs, labs)

  if (is.null(ids))
    return(kin.matrix)

  ids = as.character(ids)

  if(length(ids) != 2)
    stop2("`ids` must be a vector of length 2")
  if(!all(ids %in% labs))
    stop2("Unknown ID label: ", setdiff(ids, labs))

  kin.matrix[ids[1], ids[2]]
}


#' @rdname external_coefs
#' @export
kinship2_inbreeding = function(x, Xchrom = FALSE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  kin.matrix = kinship2_kinship(x, Xchrom = Xchrom)

  inb = numeric(pedsize(x))
  names(inb) = labels(x)

  nonf_int = nonfounders(x, internal = TRUE)
  inb[nonf_int] = kin.matrix[cbind(x$FIDX, x$MIDX)] # founders -> kin[0,0] -> excluded

  # X: males are always 1
  if(Xchrom)
    inb[getSex(x) == 1L] = 1

  inb
}
