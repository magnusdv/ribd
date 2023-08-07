#' Colourised IBD plot
#'
#' This is a pedagogical tool for illustrating the concept of
#' identity-by-descent, by representing the alleles in a pedigree by coloured
#' points or letters. By default, the alleles are placed below each pedigree
#' symbol, but this may be modified. (See examples.)
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
#' @param sep The separation between alleles within a pair, given as a multiple
#'   of the width of a pedigree symbol. Default is 1 when `pos > 0` and 0.5 for
#'   `pos = 0`.
#' @param dist The distance between pedigree symbols and the alleles, given as a
#'   multiple of symbol size. Default: 1. Ignored when `pos = 0`.
#' @param labs A logical indicating if labels should be included.
#' @param checkFounders A logical. If TRUE (default), a warning is issued if a
#'   founder has two equal alleles other than 0.
#' @param checkParents A logical. If TRUE (default), a warning is issued if
#'   someone's alleles don't match those of the parents. This a superficial test
#'   and does not catch all Mendelian errors.
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
#' # Other placements
#' ibdDraw(x, als, cols = cols, pos = c(2, 3, 1, 4))
#'
#' # Letters instead of points
#' ibdDraw(x, als, cols = cols, symbol = "text")
#'
#' # Further arguments (note that `col` is an argument of `ped.plot()`)
#' ibdDraw(x, als, cols = cols, pos = 0, symbolsize = 2.5,
#'         labs = TRUE, fill = "lightgray")
#'
#' # Mutations are warned about (unless `checkParents = FALSE`)
#' ibdDraw(x, alleles = list(1:2, 3:4, 5, 6))
#'
#'
#' ##############################
#' # Example 2: Cousin pedigree #
#' ##############################
#'
#' x = cousinPed(1) |> swapSex(3) |> relabel()
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
#' als = list(1, 2, 3, 2, 2)
#' ibdDraw(x, als, cols = 2:4)
#'
#'
#' # Restore graphics parameters
#' par(op)
#'
#' @importFrom graphics points text strheight strwidth
#' @export
ibdDraw = function(x, alleles, symbol = c("point", "text"), pos = 1, cols = NULL,
                   cex = NA, sep = 1, dist = 1, labs = FALSE, checkFounders = TRUE,
                   checkParents = TRUE, ...) {

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

  # Position
  if(!all(pos %in% 0:4))
    stop2("Illegal `pos` value: ", setdiff(pos, 0:4))
  pos = rep(pos, length.out = nInd)

  # ID labels?
  if(isTRUE(labs))
    labs = labels(x)
  else if(isFALSE(labs))
    labs = NULL

  # Alignment
  align = .pedAlignment(x, ...)

  # Add extra space to fit alleles
  addSpace = c(0,0,0,0)
  boxapprox = strwidth("ABC", "inches", cex = 1) * 2.5/3
  alleleH = strheight("M", "inches", cex = switch(symbol, point = cex*2/3, text = cex))

  botid = align$plotord[align$yall == align$yrange[2]]
  if(any(pos[botid] == 1))
    addSpace[1] = dist*boxapprox + alleleH/2

  topid = align$plotord[align$yall == align$yrange[1]]
  if(any(pos[topid] == 3))
    addSpace[3] = dist*boxapprox + alleleH/2

  leftid = align$plotord[align$xall == align$xrange[1]]
  if(any(pos[leftid] == 2))
    addSpace[2] = boxapprox*(dist+sep) + alleleH/2
  else if(any(pos[leftid] %in% c(1,3)))
    addSpace[2] = max(0, alleleH/2 - (1-sep)*boxapprox/2)

  rightid = align$plotord[align$xall == align$xrange[2]]
  print(pos[rightid])
  if(any(pos[rightid] == 4))
    addSpace[4] = boxapprox*(dist+sep) + alleleH/2
  else if(any(pos[rightid] %in% c(1,3)))
    addSpace[4] = max(0, alleleH/2 - (1-sep)*boxapprox/2)

  # Underlying pedigree plot
  p = plot(x, labs = labs, keep.par = TRUE, addSpace = addSpace, ...)

  # Height/width of ped symbols
  symh = p$scaling$boxh
  symw = p$scaling$boxw

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

    SEP = if(pos[i] > 0) sep*symw else sep*symw/2

    # Centre of allele pair
    if(pos[i] == 0) {
      X = p$x[i]
      Y = p$y[i] + symh/2
    }
    if(pos[i] == 1) {
      X = p$x[i]
      Y = p$y[i] + symh + dist*symh
    }
    else if(pos[i] == 2) {
      X = p$x[i] - symw/2 - dist*symw - SEP/2
      Y = p$y[i] + symh/2
    }
    else if(pos[i] == 3) {
      X = p$x[i]
      Y = p$y[i] - dist*symh
    }
    else if(pos[i] == 4) {
      X = p$x[i] + symw/2 + dist*symw + SEP/2
      Y = p$y[i] + symh/2
    }

    # Coordinates of pair
    if(length(als) == 2) {
      X = c(X - SEP/2, X + SEP/2)
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
