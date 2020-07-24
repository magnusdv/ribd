#' IBD triangle plot
#'
#' The IBD triangle is typically used to visualize the pairwise relatedness of
#' non-inbred individuals. Various annotations are available, including points
#' marking the most common relationships, contour lines for the kinship
#' coefficients, and shading of the unattainable region.
#'
#' For any pair of non-inbred individuals A and B, their genetic relationship
#' can be summarized by the IBD coefficients \eqn{(\kappa_0, \kappa_1,
#' \kappa_2)}{(\kappa0, \kappa1, \kappa2)}, where \eqn{\kappa_i} = P(A and B
#' share i alleles IBD at random autosomal locus). Since \eqn{\kappa_0 +
#' \kappa_1 + \kappa_2 = 1}{\kappa0 + \kappa1 + \kappa2 = 1}, any relationship
#' corresponds to a point in the triangle in the \eqn{(\kappa_0,
#' \kappa_2)}{(\kappa0, \kappa2)}-plane defined by \eqn{\kappa_0 \ge 0, \kappa_2
#' \ge 0, \kappa_0 + \kappa_2 \le 1}{\kappa0 \ge 0, \kappa2 \ge 0, \kappa0 +
#' \kappa2 \le 1}. The choice of \eqn{\kappa_0}{\kappa0} and
#' \eqn{\kappa_2}{\kappa2} as the axis variables is done for reasons of symmetry
#' and is not significant (other authors have used different views of the
#' triangle).
#'
#' As shown by Thompson (1976), points in the subset of the triangle defined by
#' \eqn{4\kappa_0\kappa_2 > \kappa_1^2}{4*\kappa0*\kappa2 > \kappa1^2} are
#' unattainable for pairwise relationships. By default this region in shaded in
#' a 'light grey' colour, but this can be modified with the `shading` argument.
#'
#' The IBD coefficients are linearly related to the kinship coefficient
#' \eqn{\phi} by the formula \deqn{\phi = 0.25\kappa_1 + 0.5\kappa_2.}{\phi =
#' 0.25*\kappa1 + 0.5*\kappa2.} By indicating values for \eqn{\phi} in the
#' `kinshipLines` argument, the corresponding contour lines are shown as dashed
#' lines in the triangle plot.
#'
#' The following abbreviations are valid entries in the `relationships`
#' argument:
#'
#' * UN = unrelated
#'
#' * PO = parent/offspring
#'
#' * MZ = monozygotic twins
#'
#' * S = full siblings
#'
#' * H,U,G = half sibling/avuncular (\strong{u}ncle)/grandparent
#'
#' * FC = first cousins
#'
#' * SC = second cousins
#'
#' * DFC = double first cousins
#'
#' * Q = quadruple first half cousins
#'
#' @param relationships A character vector indicating relationships points to be
#'   included in the plot. See Details for a list of valid entries.
#' @param kinshipLines A numeric vector (see Details).
#' @param shading The shading colour for the unattainable region.
#' @param pch Symbol used for the relationship points (see [par()]).
#' @param cexPoint A number controlling the symbol size for the relationship
#'   points.
#' @param cexText A number controlling the font size for the relationship
#'   labels.
#' @param axes A logical: Draw surrounding axis box? Default: `FALSE`.
#' @param xlab,ylab Axis labels.
#' @param cexLab A number controlling the font size for the axis labels.
#' @param xlim,ylim,mar,xpd Graphical parameters; see [par()].
#' @param keep.par A logical. If TRUE, the graphical parameters are not reset
#'   after plotting, which may be useful for adding additional annotation.
#' @return None
#' @author Magnus Dehli Vigeland
#' @references
#'
#' * E. A. Thompson (1975). _The estimation of pairwise relationships._ Annals
#' of Human Genetics 39.
#'
#' * E. A. Thompson (1976). _A restriction on the space of genetic
#' relationships._ Annals of Human Genetics 40.
#'
#' @examples
#' opar = par(no.readonly = TRUE) # store graphical parameters
#'
#' ibdTriangle()
#' ibdTriangle(kinshipLines = c(0.25, 0.125), shading = NULL, cexText = 0.8)
#'
#' par(opar) # reset graphical parameters
#'
#' @importFrom graphics abline grconvertX grconvertY layout legend mtext par
#'   plot points polygon rect segments text
#'
#' @export
ibdTriangle = function(relationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC"),
                       pch = 16, cexPoint = 1.2, cexText = 1.2,
                       kinshipLines = numeric(), shading = "lightgray",
                       xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
                       xlab = expression(kappa[0]), ylab = expression(kappa[2]),
                       cexLab = cexText, mar = c(3.1, 3.1, 1, 1), xpd = TRUE,
                       keep.par = TRUE) {

    opar = par(xpd = xpd, mar = mar, pty = "s")
    if (!keep.par)
      on.exit(par(opar))

    plot(NULL, xlim = xlim, ylim = ylim, axes = axes, ann = FALSE)

    # Axis labels
    mtext(text = c(xlab, ylab), side = 1:2, line = c(1, 0.5), las = 1, cex = cexLab)

    # impossible region shading(do borders afterwards)
    kk0 = seq(0, 1, length = 501)
    kk2 = 1 + kk0 - 2 * sqrt(kk0)
    polygon(kk0, kk2, col = shading, border = NA)
    # text(.4, .4, 'impossible region', srt = -45)

    # impossible border
    points(kk0, kk2, type = "l", lty = 3)

    # axes
    segments(c(0, 0, 0), c(0, 0, 1), c(1, 0, 1), c(0, 1, 0))

    # kinship lines
    for (phi in kinshipLines) {
        if (phi < 0 || phi > 0.5)
            stop2("kinship coefficient not in intervall [0, 0.5]", phi)
        abline(a = (4 * phi - 1), b = 1, lty = 2)
        labpos.x = 0.5 * (1.2 - (4 * phi - 1))
        labpos.y = 1.2 - labpos.x
        lab = substitute(paste(phi1, " = ", a), list(a = phi))
        text(labpos.x, labpos.y, labels = lab, pos = 3, srt = 45)
    }

    # relationships
    RELS = data.frame(
      label = c("UN", "PO", "MZ", "S", "H,U,G", "FC", "SC", "DFC", "Q"),
      k0 = c(1, 0, 0, 1/4, 1/2, 3/4, 15/16, 9/16, 17/32),
      k1 = c(0, 1, 0, 1/2, 1/2, 1/4, 1/16, 6/16, 14/32),
      k2 = c(0, 0, 1, 1/4, 0, 0, 0, 1/16, 1/32),
      pos = c(1, 1, 4, 4, 1, 1, 1, 3, 4))

    if (length(relationships) > 0) {
        rels = RELS[RELS$label %in% relationships, ]
        points(rels$k0, rels$k2, pch = pch, cex = cexPoint)
        text(rels$k0, rels$k2, labels = rels$label, pos = rels$pos, cex = cexText, offset = 0.7)
    }
}



#' Add points to the IBD triangle
#'
#' Utility function for plotting points in the IBD triangle.
#'
#' @param kappa Coordinates of points to be plotted in the IBD triangle. Valid
#'   input types are:
#'
#'   * A numerical vector of length 2 or 3. In the latter case `kappa[c(1, 3)]`
#'   is used.
#'
#'   * A matrix of data frame, whose column names must include either `k0` and
#'   `k2`, `kappa0` and `kappa2`, or `ibd0` and `ibd2`.
#'
#'   * A list (and not a data frame), in which case an attempt is made to bind
#'   the elements row-wise.
#'
#' @param new A logical indicating if a new triangle should be drawn.
#' @param col,cex,pch,lwd Parameters passed onto [points()].
#' @param labels A character of same length as the number of points, or a single
#'   logical `TRUE` or `FALSE`. If `TRUE`, an attempt is made to create labels
#'   by pasting columns `ID1` and `ID2` in `kappa`, if these exist. By default,
#'   no labels are plotted.
#' @param colLab,cexLab,pos,adj Parameters passed onto [text()] (if `labels` is
#'   non-NULL).
#' @param keep.par A logical. If TRUE, the graphical parameters are not reset
#'   after plotting, which may be useful for adding additional annotation.
#' @param \dots Plot arguments passed on to `ibdTriangle()`.
#'
#' @return None
#' @author Magnus Dehli Vigeland
#'
#' @examples
#' showInTriangle(c(3/8, 1/8), label = "3/4 siblings", pos = 1)
#'
#' @export
showInTriangle = function(kappa, new = TRUE, col = "blue",
                          cex = 1, pch = 4, lwd = 2, labels = FALSE,
                          colLab = col, cexLab = 0.8,
                          pos = 1, adj = NULL, keep.par = TRUE, ...) {

  if(is.vector(kappa) && !is.list(kappa)) {
    if(!is.numeric(kappa))
      stop2("Vector `kappa` is not numeric: ", kappa)
    if(!length(kappa) %in% 2:3)
      stop2("Vector `kappa` must have length 2 or 3. Received length: ", length(kappa))
    if(is.null(names(kappa)))
      names(kappa) = if(length(kappa) == 2) paste0("kappa", c(0, 2)) else paste0("kappa", 0:2)

    kappa = as.data.frame(as.list(kappa))
  }
  else if(is.matrix(kappa))
    kappa = as.data.frame(kappa)
  else if(is.list(kappa) && !is.data.frame(kappa))
    kappa = do.call(rbind, kappa)
  else if(!is.data.frame(kappa))
    stop2("Wrong input format of `kappa`: ", class(kappa))

  ### `kappa` is now a data frame
  df = kappa

  # Extract k0 and k2
  allowedColnames = list(c("kappa0", "kappa2"),
                         c("k0", "k2"),
                         c("ibd0", "ibd2"))
  nms = names(kappa)
  k0 = k2 = NULL
  for (colnms in allowedColnames) {
    if(all(colnms %in% nms)) {
      k0 = df[[colnms[1]]]
      k2 = df[[colnms[2]]]
      break
    }
  }
  if(is.null(k0))
    stop2("Columns names not recognised.")

  # Labels
  if(isTRUE(labels)) {
    id1 = if("ID1" %in% nms) df$ID1 else if ("id1" %in% nms) df$id1 else NULL
    id2 = if("ID2" %in% nms) df$ID2 else if ("id2" %in% nms) df$id2 else NULL
    if(is.null(id1) || is.null(id2))
      stop2("Trying to create labels, but cannot find column `id1` and `id2`")
    labels = paste(id1, id2, sep = "-")
  }


  if(is.character(labels) && length(labels) != length(k0))
    stop2("When `labels` is a character, it must have the same length as `k0`")

  if(new) {
    if(!keep.par) {
      opar = par(no.readonly = TRUE) # store graphical parameters
      on.exit(par(opar))
    }

    ibdTriangle(...)
  }

  points(k0, k2, col = col, pch = pch, lwd = lwd, cex = cex)

  if(is.character(labels))
    text(k0, k2, labels = labels, col = colLab, cex = cexLab, pos = pos, adj = adj)
}
