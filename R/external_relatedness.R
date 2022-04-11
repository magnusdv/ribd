#' Relatedness coefficients by other programs
#'
#' Wrappers for functions in other packages or external programs, computing
#' various relatedness coefficients.
#'
#' `kinship2_kinship()` and `kinship2_inbreeding()` are wrappers of
#' [kinship2::kinship()].
#'
#' `idcoefs()` wraps `identity::identity.coefs()`, which is an R interface for
#' the C program `IdCoefs` written by Mark Abney (2009). The `identity.coefs()`
#' function sometimes causes R to crash, motivating an alternative wrapper
#' `idcoefs2()` which executes an external call to the original C program
#' `IdCoefs` (version 2.1.1). For this to work, `IdCoefs` must be installed on
#' the computer (see link in the References section below) and available on the
#' system's PATH variable. The function `idcoefs2()` then writes the necessary
#' files to disk and calls `IdCoefs` via [system()].
#'
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#' @param ids A integer vector of length 2.
#' @param Xchrom A logical, indicating if the autosomal (default) or
#'   X-chromosomal coefficients should be computed.
#' @param verbose A logical, indicating if messages from `IdCoefs` should be
#'   printed.
#' @param cleanup A logical: If TRUE, the pedfile and sample file created for
#'   the `IdCoefs` run are deleted automatically.
#' @return For `kinship2_inbreeding()`, a numerical vector with inbreeding
#'   coefficients, named with ID labels.
#'
#'   For `kinship2_kinship()`, either a single numeric (if `ids` is a pair of
#'   pedigree members) or the whole kinship matrix, with the ID labels as
#'   dimnames.
#'
#'   For `idcoefs()` and `idcoefs2()`, a numerical vector of length 9 (in the
#'   standard order of Jacquard's condensed identity coefficients).
#'
#' @author Magnus Dehli Vigeland
#'
#' @seealso [kinship2::kinship()], `identity::identity.coefs()`
#'
#' @references Abney, Mark (2009). _A graphical algorithm for fast computation
#'   of identity coefficients and generalized kinship coefficients._
#'   Bioinformatics, 25, 1561-1563.
#'   <https://home.uchicago.edu/~abney/abney_web/Software.html>
#'
#' @examples
#' # A random pedigree with 2 founders and 5 matings
#' p = randomPed(g = 5, founders = 2, seed = 123)
#'
#' ### Kinship matrix
#'
#' # Autosomal: Check that ribd agrees with kinship2
#' stopifnot(identical(
#'   kinship(p),          # ribd
#'   kinship2_kinship(p)  # kinship2
#'
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
#'
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


idcoefs = function(x, ids) {
  if(!requireNamespace("identity", quietly = TRUE))
    stop2("Package `identity` must be installed for this function to work")
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must be a vector of length 2")

  ids_int = internalID(x, ids)
  ped = cbind(1:pedsize(x), x$FIDX, x$MIDX)
  j = identity::identity.coefs(ids_int, ped)[2, 3:11]

  # Swap symmetric coefficients if needed
  if(ids_int[1] > ids_int[2])
    j[c(3,4,5,6)] = j[c(5,6,3,4)]

  j
}


#' @rdname external_coefs
#' @importFrom utils read.table write.table
#' @export
idcoefs2 = function(x, ids, verbose = FALSE, cleanup = TRUE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must be a vector of length 2")
  if(!nzchar(pth <- Sys.which("idcoefs"))) stop2("Executable `idcoefs` not found")

  x = parentsBeforeChildren(x)
  ped = as.data.frame(x)[, 1:3]

  write.table(ped, file = "__paramlink2idcoefs__.ped", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(ids, file = "__paramlink2idcoefs__.sample", quote = FALSE, row.names = FALSE, col.names = FALSE)
  command = "idcoefs -p __paramlink2idcoefs__.ped -s __paramlink2idcoefs__.sample -o __paramlink2idcoefs__.output"
  run = suppressWarnings(system(command, intern = TRUE))

  if(verbose)
    message(run)
  res = read.table("__paramlink2idcoefs__.output", as.is = TRUE)

  if (cleanup)
    unlink(dir(pattern = "__paramlink2idcoefs__"))

  as.numeric(res[2, 3:11])
}
