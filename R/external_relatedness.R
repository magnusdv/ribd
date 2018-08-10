#' Relatedness coefficients by other programs
#'
#' Wrappers for functions in other packages or external programs, computing
#' various relatedness coefficients.
#'
#' Both `kinship2_inbreeding` and `kinship2_kinship` are wrappers of
#' [kinship2::kinship()]. Similarly, `jacquard` wraps
#' [identity::identity.coefs()], which is an R interface for the C program
#' `IdCoefs` written by Mark Abney (2009). The `identity.coefs()` function
#' sometimes causes R to crash, hence I have provided an alternative wrapper,
#' `jacquard2`, which executes an external call to the original C program
#' `IdCoefs` (version 2.1.1). For this to work, `IdCoefs` must be installed on
#' the computer (see link in the References section below) and the executable
#' placed in a folder included in the PATH variable. The `jacquard2` wrapper
#' works by writing the necessary files to disk and calling `IdCoefs` via
#' [system()].
#'
#' @param x A pedigree, in the form of a [`pedtools::ped`] object.
#' @param ids A integer vector of length 2.
#' @param verbose A logical, indicating if messages from IdCoefs should be
#'   printed.
#' @param cleanup A logical: If TRUE, the pedfile and sample file created for
#'   the IdCoefs run are deleted automatically.
#' @return For `kinship2_inbreeding`, a named numerical vector with the
#'   inbreeding coefficients and ID labels as names.
#'
#'   For `kinship2_kinship`, either a single numeric (if `ids` is a pair of
#'   pedigree members) or the whole kinship matrix, with the ID labels as
#'   dimnames.
#'
#'   For `jaquard` and `jaquard2`, a numerical vector of length 9 (in the
#'   standard order of Jacquard's condensed identity coefficients).
#' @author Magnus Dehli Vigeland
#' @seealso [kinship2::kinship()], [identity::identity.coefs()]
#' @references Abney, Mark (2009). _A graphical algorithm for fast computation
#'   of identity coefficients and generalized kinship coefficients._
#'   Bioinformatics, 25, 1561-1563.
#'   <http://home.uchicago.edu/~abney/abney_web/Software.html>
#'
#' @examples
#' library(pedtools)
#'
#' # Offspring of first cousins
#' x = cousinsPed(1, child=TRUE)
#' inb = kinship2_inbreeding(x)
#' stopifnot(inb[9] == 1/16)
#'
#' # if ID labels are not 1:9, care must be taken in extracting correct elements.
#' y = relabel(x, 9:1)
#' y
#' inb = kinship2_inbreeding(y)
#'
#' # The child's Id is now '1'.
#' child = leaves(y)
#'
#' # Note the difference.
#' inb[1]   #wrong
#' inb['1'] #correct
#'
#' # The inbreeding coeff of the child equals the kinship coeff of parents
#' par = parents(y, child)
#' kin = kinship2_kinship(y, par)
#' stopifnot(inb[child]==kin)
#'
#' @name relatednessCoeff
NULL

#' @rdname relatednessCoeff
#' @importFrom kinship2 kinship
#' @export
kinship2_inbreeding = function(x) {
  kin.matrix = kinship2_kinship(x)

  inb = numeric(pedsize(x))
  names(inb) = x$LABELS

  nonf_int = nonfounders(x, internal = TRUE)
  inb[nonf_int] = kin.matrix[cbind(x$FID, x$MID)] # founders -> kin[0,0] -> excluded

  inb
}

#' @rdname relatednessCoeff
#' @export
kinship2_kinship = function(x, ids = NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  if (!is.null(ids)) {
    if(length(ids) != 2) stop2("`ids` must be a vector of length 2")
    if(!all(ids %in% x$LABELS)) stop2("Unknown ID label: ", setdiff(ids, x$LABELS))
  }

  kin.matrix = kinship2::kinship(id = x$ID, dadid = x$FID, momid = x$MID)
  dimnames(kin.matrix) = list(x$LABELS, x$LABELS)

  if (is.null(ids))
    return(kin.matrix)
  kin.matrix[as.character(ids[1]), as.character(ids[2])]
}

#' @rdname relatednessCoeff
#' @export
jacquard = function(x, ids) {
  if (!requireNamespace("identity", quietly = TRUE))
      stop2("Package `identity` must be installed for this function to work")
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must be a vector of length 2")

  ids_int = internalID(x, ids)
  ped = cbind(x$ID, x$FID, x$MID)
  j = identity::identity.coefs(ids_int, ped)[2, 3:11]

  # Swap symmetric coefficients if needed
  if(ids_int[1] > ids_int[2])
    j[c(3,4,5,6)] = j[c(5,6,3,4)]

  j
}


#' @rdname relatednessCoeff
#' @importFrom utils read.table write.table
#' @export
jacquard2 = function(x, ids, verbose = FALSE, cleanup = TRUE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(ids) != 2) stop2("`ids` must be a vector of length 2")

  x = parents_before_children(x)
  ped = as.data.frame(x)[, 1:3]

  write.table(ped, file = "__paramlink2idcoefs__.ped", quote = F, row.names = F, col.names = F)
  write.table(ids, file = "__paramlink2idcoefs__.sample", quote = F, row.names = F, col.names = F)
  command = "idcoefs -p __paramlink2idcoefs__.ped -s __paramlink2idcoefs__.sample -o __paramlink2idcoefs__.output"
  run = suppressWarnings(system(command, intern = T))

  if (verbose)
    print(run)
  res = read.table("__paramlink2idcoefs__.output", as.is = T)

  if (cleanup)
    unlink(dir(pattern = "__paramlink2idcoefs__"))

  as.numeric(res[2, 3:11])
}
