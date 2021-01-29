#' Colourised IBD plot
#'
#' This is a pedagogical tools for illustrating the concept of
#' identity-by-descent, by representing the alleles in a pedigree by coloured
#' points or letters. By default, the alleles are placed below each pedigree
#' symbols, but any positions are possible, including inside. (See examples.)
#'
#' @param x A `ped` object.
#' @param alleles A list of length `pedsize(x)`. Each element should consist of
#'   one or two integers, representing different colours. Zeroes produce
#'   "greyed-out" alleles.
#' @param symbol Either "point" or "text".
#' @param pos A vector recycled to the length of `labels(x)`, indicating allele
#'   placement relative to the pedigree symbols: 0 = inside; 1 = below; 2 =
#'   left; 3 = above; 4 = right. By default, all are placed below.
#' @param cols A colour vector corresponding to the integers occurring in
#'   `alleles`.
#' @param cex An expansion factor for the allele points/letters. Default: 3 for
#'   points and 2 for text.
#' @param sep The separation between haplotypes within a pair, given as a
#'   multiple of the width of a pedigree symbol. Default: 0.5 when `pos = 0` and
#'   1 otherwise.
#' @param dist The distance between pedigree symbols and the alleles, given as a
#'   multiple of the height of a pedigree symbol. Default: 1. Ignored when `pos
#'   = 0`.
#' @param labs A logical indicating if labels should be included.
#' @param checkFounders A logical. If TRUE (default), a warning is issued if a
#'   founder has two equal alleles other than 0.
#' @param checkParents A logical. If TRUE (default), a warning is issued if
#'   someone's alleles don't match those of the parents. This a superficial test
#'   and does not catch all Mendelian errors.
#' @param margin Plot margins (bottom, left, top, right).
#' @param ... Further arguments passed on to `plot.ped()`.
#'
#' @return The plot structure is returned invisibly.
#'
#' @seealso [pedtools::plot.ped()], `ibdsim2::haploDraw()`
#'
#' @examples
#' op = par(no.readonly = TRUE)
#'
#' ###############################
#' # Example 1: A family quartet #
#' ###############################
#'
#' x = nuclearPed(2)
#' als = list(1:2, 3:4, c(1,3), c(2,3))
#'
#' # Default options
#' ibdDraw(x, als)
#'
#' # Nicer colors
#' cols = c(7, 3, 2, 4)
#' ibdDraw(x, als, cols = cols)
#'
#' # Inside the pedigree symbols
#' ibdDraw(x, als, cols = cols, pos = 0, symbolsize = 2.5)
#'
#' # Other placements (margins depend on device - may need adjustment)
#' ibdDraw(x, als, cols = cols, pos = c(2, 4, 1, 1),
#'         margin = c(2, 6, 2, 6))
#'
#' # Letters instead of points
#' ibdDraw(x, als, cols = cols, symbol = "text", cex = 2)
#'
#' # Further arguments (note that `col` is an argument of `ped.plot()`)
#' ibdDraw(x, als, cols = cols, pos = 0, symbolsize = 2,
#'         labs = TRUE, hatched = 3:4, col = "blue")
#'
#' # Mutations are warned about (unless `checkParents = FALSE`)
#' ibdDraw(x, alleles = list(1:2, 3:4, 5, 6))
#'
#'
#' ##############################
#' # Example 2: Cousin pedigree #
#' ##############################
#'
#' x = swapSex(cousinPed(1), 3)
#' als = list(1:2, 3:4, NULL, c(1,3), c(2,3), NULL, 3, 3)
#'
#' cols = c(7, 3, 2, 4)
#' ibdDraw(x, als, cols = cols, dist = 0.8)
#' ibdDraw(x, als, cols = cols, dist = 0.8, symbol = "text")
#'
#' # Alternative: 0's give greyed-out alleles
#' als2 = list(1:2, 3:4, c(0,0), c(1,3), c(2,3), c(0,0), c(0,3), c(3,0))
#'
#' ibdDraw(x, als2, cols = cols, dist = 0.8)
#' ibdDraw(x, als2, cols = cols, dist = 0.8, symbol = "text")
#'
#'
#' ############################
#' # Example 3: X inheritance #
#' ############################
#'
#' x = nuclearPed(2, sex = c(1, 2))
#' als = list(1, 2:3, 3, c(1, 3))
#' ibdDraw(x, als, cols = c(3, 7, 2))
#'
#'
#' #################################
#' # Example 4: mtDNA inheritance  #
#' #################################
#'
#' x = linearPed(2, sex = 2)
#' als = list(1, 2, 2, 3, 2)
#' ibdDraw(x, als, cols = 2:4)
#'
#'
#' # Restore graphics parameters
#' par(op)
#'
#' @importFrom graphics points text
#' @export
ibdDraw = function(x, alleles, symbol = c("point", "text"), pos = 1, cols = NULL,
                   cex = NA, sep = NULL, dist = 1, labs = FALSE, checkFounders = TRUE,
                   checkParents = TRUE, margin = c(1, 1, 1, 1), ...) {

  if(!is.ped(x))
    stop2("Argument `x` must be a `ped` object")

  nInd = pedsize(x)

  # Fix `alleles`
  if(!is.null(names(alleles))) {
    a = vector(nInd, mode = "list")
    a[internalID(x, names(alleles))] = alleles
    alleles = a
  }

  if(length(alleles) < nInd)
    alleles = c(alleles, vector(nInd - length(alleles), mode = "list"))

  uniqAls = unique.default(unlist(alleles))

  if(!setequal(setdiff(uniqAls, 0), 1:max(uniqAls)))
    stop2("Consecutive numbers should be used as alleles. Missing: ",
          setdiff(1:max(uniqAls), uniqAls))


  symbol = match.arg(symbol)

  if(is.na(cex))
    cex = switch(symbol, point = 3, text = 2)

  # Default colours
  if(is.null(cols))
    cols = 1:max(uniqAls)
  else if(length(cols) < max(uniqAls))
    stop2("Too few colours")

  # Position and separation
  pos = rep(pos, length.out = nInd)
  if(is.null(sep))
    sep = ifelse(pos == 0, 0.5, 1)
  else
    sep = rep(sep, length.out = nInd)

  # ID labels?
  if(isTRUE(labs))
    labs = labels(x)
  else if(isFALSE(labs))
    labs = NULL

  # Underlying pedigree plot
  p = plot(x, labs = labs, keep.par = TRUE, margin = margin, ...)

  # Height/width of ped symbols
  symh = p$boxh
  symw = p$boxw
  SEP = sep * symw
  DIST = dist * symh

  fou = founders(x, internal = TRUE)

  # Loop through all individuals in pedigree
  for(i in 1:nInd) {
    als = alleles[[i]]
    if(is.null(als))
      next

    # Various checks
    lab = labels(x)[i]

    if(length(als) > 2)
      stop2("Individual has more than two alleles: ", lab)

    if(checkFounders && i %in% fou) {
      if(length(als) == 2 && all(als > 0) && als[1] == als[2])
        warning("Founder ", lab, " has two identical alleles", call. = FALSE)
    }

    if(checkParents && !i %in% fou) {
      fa = alleles[[father(x, i, internal = TRUE)]]
      mo = alleles[[mother(x, i, internal = TRUE)]]
      if(!is.null(fa) && !is.null(mo) && !all(als %in% c(0, fa, mo)))
        warning("Allele of individual ", labels(x)[i], " unseen in parents: ",
                setdiff(als, c(0, fa, mo)), call. = FALSE)
    }

    # Centre of allele pair
    if(pos[i] == 0) {
      X = p$x[i]
      Y = p$y[i] + symh/2
    }
    if(pos[i] == 1) {
      X = p$x[i]
      Y = p$y[i] + symh + DIST
    }
    else if(pos[i] == 2) {
      X = p$x[i] - symw/2 - DIST - SEP[i]/2
      Y = p$y[i] + symh/2
    }
    else if(pos[i] == 3) {
      X = p$x[i]
      Y = p$y[i] - DIST
    }
    else if(pos[i] == 4) {
      X = p$x[i] + symw/2 + DIST + SEP[i]/2
      Y = p$y[i] + symh/2
    }

    # Coordinates of pair
    if(length(als) == 2) {
      X = c(X - SEP[i]/2, X + SEP[i]/2)
      Y = c(Y,Y)
    }

    if(symbol == "point") {
      col = ifelse(als > 0, 1, 8)
      bg = ifelse(als > 0, cols[als], 0)
      points(X, Y, pch = 21, cex = cex, col = col, bg = bg)
    }
    else if(symbol == "text") {
      txt = ifelse(als > 0, LETTERS[als], "*")
      col = ifelse(als > 0, cols[als], 8)
      text(X, Y, txt, col = col, cex = cex)
    }
  }

  invisible(p)
}

#' ibdDraw(nuclearPed(), al = list(1:2, 3:4, c(1,3)), cols = c(3,7))
