#' IBD triangle plot
#'
#' A relationship triangle used to used to visualize the kappa coefficients of
#' between non-inbred individuals. Various annotations are available, including
#' points marking the most common relationships, contour lines for the kinship
#' coefficients, and shading of the unattainable region. The companion function
#' [showInTriangle()] (which plots user-specified points onto the triangle) is
#' probably the most useful for end users.
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
#' a light grey colour, but this can be modified with the `shading` argument.
#'
#' The IBD coefficients are linearly related to the kinship coefficient
#' \eqn{\varphi} by the formula \eqn{\varphi = 0.25\kappa_1 + 0.5\kappa_2.}{\phi =
#' 0.25*\kappa1 + 0.5*\kappa2.} By indicating values for \eqn{\varphi} in the
#' `kinshipLines` argument, the corresponding contour lines are shown in the
#' triangle plot. (Currently only when `plotType = "base"`.)
#'
#' @param relationships A character vector indicating the *fixed* relationships
#'   points to be included in the plot. Valid entries are those in the `label`
#'   column of [basicRelationships].
#' @param plotType Either "base" (default), "ggplot2" or "plotly". Abbreviations
#'   are allowed.
#' @param kinshipLines A numeric vector (see Details).
#' @param shortLines A logical indicating if the kinship lines (if present)
#'   should be restricted to the interior of the triangle.
#' @param shading The shading colour for the unattainable region.
#' @param pch Symbol used for the relationship points (see [par()]).
#' @param cexPoint A number controlling the symbol size for the relationship
#'   points.
#' @param cexText A number controlling the font size for the relationship
#'   labels.
#' @param axes A logical: Draw surrounding axis box? Default: `FALSE`.
#' @param xlab,ylab Axis labels.
#' @param title Main title (absent by default).
#' @param cexAxis A number controlling the font size for the axis labels.
#' @param cexLab Deprecated; use `cexAxis` instead.
#' @param xlim,ylim,xpd,las Graphical parameters; see [par()]. (Base plot only.)
#' @param mar Graphical parameter; see [par()]. For ggplot2, this is ignored
#'   unless it is a ggplot2::margin() object.
#' @param keep.par A logical. If TRUE, the graphical parameters are not reset
#'   after plotting, which may be useful for adding additional annotation. (Base
#'   plot only.)
#' @param \dots Further arguments; currently not used.
#'
#' @return `NULL` if `plotType = 'base'`; otherwise the plot object.
#'
#' @references
#'
#' * E. A. Thompson (1975). _The estimation of pairwise relationships._ Annals
#' of Human Genetics 39.
#'
#' * E. A. Thompson (1976). _A restriction on the space of genetic
#' relationships._ Annals of Human Genetics 40.
#'
#' @seealso [showInTriangle()], [kappaIBD()]
#'
#' @examples
#' opar = par(no.readonly = TRUE) # store graphical parameters
#'
#' ibdTriangle()
#' ibdTriangle(kinshipLines = c(0.25, 0.125), shading = NULL, cexText = 0.7)
#' ibdTriangle(kinshipLines = c(0.25, 0.125), shortLines = TRUE, pch = 15)
#' ibdTriangle(relationships = c("UN", "PO", "MZ", "S"),
#'             xlab = "k0", ylab = "k2", las = 0, axes = TRUE, cexAxis =1.6)
#'
#' par(opar) # reset graphical parameters
#'
#' @importFrom graphics abline grconvertX grconvertY layout legend mtext par
#'   plot points polygon rect segments text
#'
#' @export
ibdTriangle = function(relationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC"),
                       plotType = c("base", "ggplot2", "plotly"), pch = 19,
                       cexPoint = 1.2, cexText = 1.2, cexAxis = cexText,
                       kinshipLines = numeric(), shortLines = FALSE, shading = "gray90",
                       xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, las = 1,
                       xlab = expression(kappa[0]), ylab = expression(kappa[2]),
                       title = NULL, mar = c(2.1, 2.1, 1, 1), xpd = TRUE, keep.par = TRUE,
                       cexLab = NULL, ...) {

  if(!is.null(cexLab)) {
    message("Argument `cexLab` of `ibdTriangle()` has been renamed to `cexAxis`")
    cexAxis = cexLab
  }

  plotType = match.arg(plotType)
  if(plotType == "ggplot2") {
    return(ibdTriangleGG(relationships = relationships, pch = pch,
                         cexPoint = cexPoint, cexText = cexText,
                         cexAxis = cexAxis, shading = shading,
                         xlim = xlim, ylim = ylim, las = las,
                         xlab = xlab, ylab = ylab, title = title, mar = mar, ...))
  }
  if(plotType == "plotly") {
    return(ibdTrianglePlotly(relationships = relationships, cexPoint = cexPoint,
                         cexText = cexText, cexAxis = cexAxis, shading = col2hex(shading),
                         xlab = xlab, ylab = ylab, ...))
  }


  # Base plot ---------------------------------------------------------------

  opar = par(xpd = xpd, mar = mar, pty = "s")
  if (!keep.par)
    on.exit(par(opar), add = TRUE)

  plot(NULL, xlim = xlim, ylim = ylim, axes = axes, ann = FALSE)

  # Axis labels
  mtext(text = c(xlab, ylab), side = 1:2, line = c(1, 0.5), las = las, cex = cexAxis)

  # Impossible region shading with dotted border
  t = seq(0, 1, length = 201)
  polygon(t^2, (1-t)^2, col = shading, border = NA)
  points(t^2, (1-t)^2, type = "l", lty = 3)
  # text(.4, .4, 'impossible region', srt = -45)

  # axes
  segments(c(0, 0, 0), c(0, 0, 1), c(1, 0, 1), c(0, 1, 0))

  # kinship lines
  for (phi in kinshipLines) {
    if (phi < 0 || phi > 0.5)
      stop2("kinship coefficient not in intervall [0, 0.5]", phi)

    if(shortLines) {
      if(phi > 1/4)
        segments(x0 = 0, y0 = 4*phi-1, x1 = 1-2*phi, y1 = 2*phi, lty = 2)
      else
        segments(x0 = 1-4*phi, y0 = 0, x1 = 1-2*phi, y1 = 2*phi, lty = 2)

      labpos.x = 1-2*phi + 0.025
      labpos.y = 2*phi + 0.025
      pos = NULL; adj = c(0, 0.5)
    }
    else {
      abline(a = (4 * phi - 1), b = 1, lty = 2)
      labpos.x = 0.5 * (1.2 - (4 * phi - 1))
      labpos.y = 1.2 - labpos.x
      pos = 3; adj = NULL
    }
    lab = substitute(paste(phi1, " = ", a), list(a = phi))
    text(labpos.x, labpos.y, labels = lab, pos = pos, adj = adj, srt = 45, cex = cexText)
  }

  # Fixed relationships
  allrels = ribd::basicRelationships
  rels = allrels[allrels$label %in% relationships, , drop = FALSE]
  if(nrow(rels)) {
    points(rels$kappa0, rels$kappa2, pch = pch, cex = cexPoint)
    text(rels$kappa0, rels$kappa2, labels = rels$label, pos = rels$pos,
         cex = cexText, offset = 0.7)
  }

  invisible()
}


ibdTriangleGG = function(relationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC"),
                         pch = 19, cexPoint = 1.2, cexText = cexPoint,
                         cexAxis = cexText, shading = "gray90",
                         xlim = c(0,1), ylim = c(0,1), las = 1, mar = NULL,
                         xlab = expression(kappa[0]), ylab = expression(kappa[2]),
                         title = NULL, ...) {

  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop2("Package `ggplot2` must be installed for this to work")

  # Fixed
  allrels = ribd::basicRelationships
  rels = allrels[allrels$label %in% relationships, , drop = FALSE]

  # Tweak label pos
  adj = 0.02
  pos = rels$pos
  rels$x.txt = rels$kappa0 + ifelse(pos == 2, -2*adj, ifelse(pos == 4, 2*adj, 0))
  rels$y.txt = rels$kappa2 + ifelse(pos == 1, -adj, ifelse(pos == 3, adj, 0))
  rels$hjust = ifelse(pos == 2, 1, ifelse(pos == 4, 0, 0.5))
  rels$vjust = ifelse(pos == 1, 1, ifelse(pos == 3, 0, 0.5))

  # Impossible region shading(do borders afterwards)
  t = seq(0, 1, length = 201)
  imp = data.frame(kappa0 = t^2, kappa2 = (1 - t)^2)

  # Main triangle
  p = ggplot2::ggplot(rels, ggplot2::aes(kappa0, kappa2)) +
    ggplot2::geom_polygon(data = imp, fill = shading, color = "black", linetype = 3) +
    ggplot2::geom_polygon(data = data.frame(kappa0 = c(0, 1, 0), kappa2 = c(0, 0, 1)),
                          fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.5 * cexPoint, shape = pch) +
    ggplot2::geom_text(ggplot2::aes(x = `x.txt`, y = `y.txt`, label = label, hjust = hjust,
                                    vjust = vjust), size = 3.88 * cexText) +
    ggplot2::coord_fixed(ratio = 1, clip = "off", xlim = xlim, ylim = ylim) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.title   = ggplot2::element_text(size = 11 * cexAxis),
                   axis.title.y = ggplot2::element_text(angle = if(las == 0) 90 else 0),
                   plot.margin = if(inherits(mar, "margin")) mar else NULL,
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 14 * cexAxis))

  p
}


ibdTrianglePlotly = function(relationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC"),
                         cexPoint = 1.2, cexText = cexPoint, cexAxis = cexText,
                         shading = "gray90", xlab = "k0", ylab = "k2", ...) {

  if(!requireNamespace("plotly", quietly = TRUE))
    stop2("Package `plotly` must be installed for this to work")

  # Fixed relationships coordinates
  allrels = ribd::basicRelationships
  rels = allrels[allrels$label %in% relationships, , drop = FALSE]

  # Adjust label positions
  plotlyPos = c(UN = "bottom right", PO = "bottom center", MZ = "middle right",
                S = "middle right", H = "bottom center", A = "bottom center",
                G = "bottom center", `H,U,G` = "bottom center", FC = "bottom center",
                SC = "bottom center", DFC = "top right", Q = "middle right")

  rels$pos = plotlyPos[rels$label]
  rels$x = rels$kappa0 + ifelse(rels$pos == "middle right", 1/50, 0)
  rels$y = rels$kappa2 - ifelse(startsWith(rels$pos, "bottom"), 1/50, 0)

  # Illegal region
  t = seq(0, 1, length = 101)
  xIlleg = t^2
  yIlleg = (1-t)^2

  # Axis labels
  xlab = .fixPlotlyText(xlab)
  ylab = .fixPlotlyText(ylab)

  # Plot: static part
  p = plotly::plot_ly() |>
    plotly::config(displaylogo = FALSE, displayModeBar = FALSE) |>
    plotly::layout(
      margin = list(l = 20, r = 20, b = 20, t = 20),
      xaxis = list(range = c(-0.05, 1.1), visible = TRUE, showticklabels = FALSE,
                   zeroline = FALSE, showgrid = FALSE, showline = FALSE,
                   title = list(text = xlab, standoff = 0)),
      yaxis = list(range = c(-0.1, 1.05), visible = TRUE, showticklabels = FALSE,
                   zeroline = FALSE, showgrid = FALSE, showline = FALSE,
                   scaleanchor = "x", scaleratio = 1,
                   title = list(text = ylab, standoff = 0))) |>
    plotly::add_segments(x = c(0, 0, 0), y = c(0, 0, 1), xend = c(1, 0, 1),
                         yend = c(0, 1, 0), hoverinfo = 'skip',
                         line = list(color = "black", width = 1),
                         showlegend = FALSE) |>
    plotly::add_polygons(x = xIlleg, y = yIlleg, fillcolor = col2hex(shading),
                         line = list(color = "black", width = 0.5, dash = "dot"),
                         hoverinfo = 'skip', showlegend = FALSE) |>
    plotly::add_markers(data = rels, x = ~kappa0, y = ~kappa2, hoverinfo = 'skip',
                        marker = list(color = "black", size = 7*cexPoint),
                        showlegend = FALSE) |>
    plotly::add_text(data = rels, x = ~x, y = ~y, text = ~label, size = I(15*cexText),
                     color = I("black"), textposition = ~pos, hoverinfo = 'skip',
                     showlegend = FALSE)

  if(inherits(xlab, "TeX") || inherits(ylab, "TeX"))
    p = p |> plotly::config(mathjax = 'cdn')

  p
}

.fixPlotlyText = function(s) {
  if(is.expression(s))
    s = deparse(s[[1]])
  if(s == "kappa[0]")
    s = plotly::TeX("\\kappa_0")
  if(s == "kappa[2]")
    s = plotly::TeX("\\kappa_2")
  s
}


#' Add points to the IBD triangle
#'
#' Utility function for plotting kappa coefficients in the IBD triangle. This
#' was previously only implemented as a base R plot, but is now also available
#' in `ggplot2` and `plotly` formats, controlled by the argument `plotType`.
#' Labels are often easier to read in the two latter versions: The `ggplot2`
#' version uses `ggrepel` to separate labels, while `plotly` enables interactive
#' exploration of the plot.
#'
#' @param kappa Coordinates of points to be plotted in the IBD triangle. Valid
#'   input types are:
#'
#'   * A numerical vector of length 2 (kappa0, kappa2) or 3 (kappa0, kappa1,
#'   kappa2). In the latter case kappa1 is ignored.
#'
#'   * A matrix of data frame, whose column names must include either `k0` and
#'   `k2`, `kappa0` and `kappa2`, or `ibd0` and `ibd2`.
#'
#'   * A list (and not a data frame), in which case an attempt is made to bind
#'   the elements row-wise.
#'
#' @param plotType Either "base" (default), "ggplot2" or "plotly". Abbreviations
#'   are allowed.
#' @param new A logical indicating if a new triangle should be drawn.
#' @param ped A pedigree to be drawn in the upper right corner of the plot.
#'   Default: NULL. This only works when `plotType` is `base` or `ggplot2`.
#' @param pedBL A vector of length two, with entries in `[0,1]`, indicating the
#'   coordinates of the bottom left corner. Default: `c(0.5, 0.5)`.
#' @param pedArgs Plotting arguments for the inset pedigree. Default: NULL.
#' @param col,cex,pch,lwd Parameters controlling the appearance of points.
#' @param jitter A logical. If NULL (default), jittering is disabled for
#'   `plotType`'s "base" or "ggplot2" and enabled for "plotly".
#' @param labels A character of same length as the number of points, or a single
#'   logical `TRUE` or `FALSE`. If `TRUE`, labels are created by pasting columns
#'   `id1` and `id2` in `kappa` (if these exist) separated by `labSep`. By
#'   default, labels are disabled when `plotType = "base"`, enabled if `plotType
#'   = "ggplot2"` and enabled (interactively) if `plotType = "plotly"`.
#' @param labSep A string, by default "-".
#' @param colLab,cexLab,pos,adj Parameters controlling the appearance of labels.
#'   Ignored when `plotType = "plotly"`.
#' @param keep.par A logical. If TRUE (and `plotType = "base"`), the graphical
#'   parameters are not reset after plotting, which may be useful for adding
#'   additional annotation.
#' @param \dots Plot arguments passed on to `ibdTriangle()`.
#'
#' @return If `plotType = 'base'`, the function returns NULL; otherwise the plot
#'   object.
#'
#' @seealso [ibdTriangle()], [kappaIBD()]
#'
#' @examples
#' showInTriangle(c(3/8, 1/8), label = "3/4 siblings", pos = 1)
#'
#' # With inset pedigree
#' x = doubleCousins(1, 0, half2 = TRUE)
#' showInTriangle(c(3/8, 1/8), label = "3/4 siblings", pos = 1,
#'                ped = x, pedArgs = list(hatched = 6:7))
#'
#' # All pairs
#' k = kappaIBD(x)
#' showInTriangle(k, labels = TRUE, pos = 1:4, ped = x)
#'
#' # With jitter and variable colors
#' showInTriangle(k, labels = TRUE, pos = 1:4, jitter = TRUE, col = 1:7, ped = x)
#'
#' # Separate labels (requires ggplot2 + ggrepel)
#' # showInTriangle(k, plot = "ggplot2", col = 2:8, ped = x)
#'
#' # Interactive plot (requires plotly)
#' # showInTriangle(k, plot = "plotly", col = 2:8, pch = 0)
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats runif
#' @export
showInTriangle = function(kappa, plotType = c("base", "ggplot2", "plotly"),
                          new = TRUE, ped = NULL, pedBL = c(.5,.5), pedArgs = NULL,
                          col = 6, cex = 1, pch = 4,
                          lwd = 2, jitter = NULL, labels = NULL, colLab = col,
                          cexLab = 0.8, labSep = "-",
                          pos = 1, adj = NULL, keep.par = TRUE, ...) {

  plotType = match.arg(plotType)

  # Pedigree inset ----------------------------------------------------------

  if(!is.null(ped)) {
    if(plotType == "plotly")
      stop2("Inset pedigrees (with the `ped` argument) does not work with plotly")

    op = par(fig = c(pedBL[1], .99, pedBL[2], .99), new = FALSE, pty = "m")
    on.exit(par(op), add = TRUE)
    tryCatch(do.call(plot, args = c(list(x = ped), pedArgs)),
             error = function(e) message(paste("No pedigree drawn.\nError message:",
                                               conditionMessage(e))))

    par(fig = c(0,1,0,1), new = TRUE, pty = "s")

    p = showInTriangle(kappa, plotType, ped = NULL, new = TRUE, col = col, cex = cex,
                     pch = pch, lwd = lwd, jitter = jitter, labels = labels,
                     colLab = colLab, cexLab = cexLab, labSep = labSep,
                     pos = pos, adj = adj, keep.par = keep.par, ...)

    if(plotType == "ggplot2")
      print(p, newpage = FALSE)
    return(invisible())

  }

  # Main triangle -----------------------------------------------------------

  # Only relevant for base plot
  if(!keep.par) {
    opar = par(no.readonly = TRUE) # store graphical parameters
    on.exit(par(opar), add = TRUE)
  }

  # Draw triangle triangle
  if(new)
    p = ibdTriangle(plotType = plotType, ...)
  else if(plotType != "base")
    stop2("Argument `new = FALSE` is only meaningful when `plotType = 'base'")


  # Prepare plot data -------------------------------------------------------

  if(is.null(kappa) || (!is.null(nrow(kappa)) && nrow(kappa) == 0))
    return(switch(plotType, base = invisible(), p))

  if(is.null(labels))
    labels = switch(plotType, base = FALSE, ggplot2 = TRUE, plotly = TRUE)

  if(is.null(jitter))
    jitter = switch(plotType, base = FALSE, ggplot2 = FALSE, plotly = TRUE)

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
  df = as.data.frame(kappa)
  nms = names(df)

  # Columns with coordinates --> copy to .k0, .k2
  for (kcols in list(c("k0", "k2"), c("kappa0", "kappa2"), c("ibd0", "ibd2"))) {
    if (all(kcols %in% nms)) {
      df$.k0 = df$.k0ex = df[[kcols[1]]]
      df$.k2 = df$.k2ex = df[[kcols[2]]]
      break
    }
  }
  if(is.null(df$.k0))
    stop2("Columns names not recognised")

  # Jitter?
  if(jitter) {
    df$.k0 = df$.k0 + runif(nrow(df), -0.015, 0.015)
    df$.k2 = df$.k2 + runif(nrow(df), -0.015, 0.015)
  }

  df$col = rep_len(col2hex(col), nrow(df))
  df$pch = rep_len(pch, nrow(df))

  # Create labels from ID columns
  if (isTRUE(labels)) {
    for (idcols in list(c("id1", "id2"), c("ID1", "ID2"))) {
      if (all(idcols %in% names(df))) {
        df$.id1 = df[[idcols[1]]]
        df$.id2 = df[[idcols[2]]]
        labels = paste(df[[idcols[1]]], df[[idcols[2]]], sep = labSep)
        break
      }
    }
  }

  if(is.character(labels))
    df$.labs = labels
  else
    labels = FALSE


  # Base plot----------------------------------------------------------------

  if(plotType == "base") {

    # Points to be drawn onto the triangle
    points(df$.k0, df$.k2, col = col, pch = pch, lwd = lwd, cex = cex)

    # Labels
    if(!isFALSE(labels))
      text(df$.k0, df$.k2, labels = df$.labs, col = colLab, cex = cexLab,
           pos = pos, adj = adj)

    return(invisible())
  }


  # ggplot2 -----------------------------------------------------------------

  if (plotType == "ggplot2") {

    # Triangle
    p = p +
      ggplot2::geom_point(data = df, ggplot2::aes(.k0, .k2), color = df$col,
                          size = 2 * cex, shape = pch, stroke = sqrt(lwd))

    if(!isFALSE(labels)) {

      # Use ggrepel to add labels
      if(!requireNamespace("ggrepel", quietly = TRUE))
        stop2("Package `ggrepel` must be installed for this option to work")

      # We want the labels to repel away from each other, but also from the triangle itself.
      # The following compiles a data frame of "dummy" points tracing the static triangle data

      v = seq(0, 1, length = 40)
      traceParts = list(
        cbind(.k0 = v, .k2 = 0),                     # x axis
        cbind(.k0 = 0, .k2 = v[1:10]),               # y axis
        cbind(.k0 = 1 - v[1:10], .k2 = v[1:10]),     # diagonal
        cbind(.k0 = (1-v[1:25])^2, .k2 = v[1:25]^2), # curve
        cbind(.k0 = c(-0.03, 0, 0.03,  0.45, 0.5, 0.55,  0.72, 0.75, 0.78,
                        0.91, 0.94, 0.97, 1, 1.03) |> rep(2),
              .k2 = rep(c(-0.02, -0.04), each = 14)),  # PO, HUG, UN
        cbind(.k0 = rep(c(0.25, 0.28, 0.31), 2),
              .k2 = c(0.24, 0.26))                     # S
      )

      traceTriangle = do.call(rbind, traceParts) |> as.data.frame()

      dfRepel = rbind(df[c(".k0", ".k2", ".labs", "col")],
                      cbind(traceTriangle, .labs = "", col = ""))

      # To see the dummy repellent points:
      # p = p + ggplot2::geom_point(data = traceTriangle, ggplot2::aes(.k0, .k2), col = 3)

      p = p +
        ggrepel::geom_text_repel(ggplot2::aes(.k0, .k2, label = .labs), color = dfRepel$col,
                                 size = 3.88 * cex, data = dfRepel, max.overlaps = Inf)
    }

    return(p)
  }


  # Plotly ------------------------------------------------------------------

  if(plotType == "plotly") {

    # Add row number
    df$idx = seq_len(nrow(df))

    # Exact kappas (for labels)
    if(!isFALSE(labels))
      df$.labs = paste0("ID1: ", df$.id1, "<br>", "ID2: ", df$.id2)

    # Add interactive points
    p = p |>
      plotly::add_markers(data = df, x = ~.k0, y = ~.k2, customdata = ~idx,
                          marker = list(symbol = ~pch, color = ~col, size = ~10*cex,
                                        line = list(color = 1, width = 1)),
                          text= ~ .labs,
                          hoverinfo = "text",
                          showlegend = FALSE)
    p
  }

}

col2hex = function(col) {
    m = col2rgb(col)
    rgb(red = m[1,], green = m[2,], blue = m[3,], maxColorValue = 255)
}

utils::globalVariables(c("kappa0", "kappa2", "x.txt", "y.txt", "label", "hjust", "vjust", ".k0", ".k2", ".labs"))

