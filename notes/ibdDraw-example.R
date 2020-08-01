# Include as example in ibdDraw() when new `pedtools` is on CRAN

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
