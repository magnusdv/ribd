#' Generalised kinship pattern
#'
#' @param x A `ped` object
#' @param pattern A list of vectors of ID labels.
#' @param distinct A logical indicating if different blocks always have distinct
#'   alleles
#'
#' @return An object of class `kinpat`.
#'
#' @examples
#' kinpat(nuclearPed(2), list(1, 3:4))
#'
#' @export
kinpat = function(x, pattern, distinct = TRUE) {
  if(!is.ped(x))
    stop2("First argument must be a `ped` object")

  if(!is.list(pattern))
    stop2("`pattern` must be a list")

  names(pattern) = NULL # to avoid problems when unlisting
  ids = unlist(pattern)
  idsInt = internalID(x, ids)

  # Deterministic?
  nms = names(ids)
  detailed = !is.null(nms)
  if(detailed) {
    nmsNum = match(nms, c("", "p", "m")) - 1L  # 0, 1 and 2
    if(anyNA(nmsNum))
      stop2("Names (to indicate deterministic sampling) must be 'p' or 'm': ",
            setdiff(nms, c("", "p", "m")))
  }

  # Group index
  grp = rep(sq <- seq_along(pattern), lengths(pattern))

  # Create kinpat object
  pat = lapply(sq, function(i) {
    idx = grp == i
    g = idsInt[idx]
    if(detailed)
      g = 10L*g + nmsNum[idx]
    g
  })

  structure(pat, labs = labels(x), detailed = detailed, distinct = distinct, class = "kinpat")
}


kinpat2string = function(kp, detailed = isDetailed(kp)) {

  labs = attr(kp, "labs")

  grps = vapply(kp, function(g) {
    if(detailed) {
      lb = labs[g %/% 10]
      par = g %% 10
      lb[par > 0] = paste(lb[par > 0], c("p", "m")[par], sep = ":")
    }
    else {
      lb = labs[g]
    }
    sprintf("(%s)", paste0(lb, collapse = ","))
  }, FUN.VALUE = "")

  paste0(grps, collapse = if(attr(kp, "distinct")) " | " else " / ")
}

#' @export
print.kinpat = function(x, ...) {
  cat(kinpat2string(x), "\n", sep = "")
}

kinpat2list = function(kp) {
  labs = attr(kp, "labs")
  detailed = isDetailed(kp)

  lapply(kp, function(g) {
    if(detailed) {
      s = labs[g %/% 10]
      names(s) = c("", "p", "m")[g %% 10 + 1]
      s
    }
    else
      labs[g]
  })
}

kinpatSort = function(kp) {

  # Function for sorting a single group
  sortGroup = function(g) .mysort(g, decreasing = TRUE)

  # Quick return if just one group
  if(length(kp) == 1) {
    kp[[1]] = sortGroup(kp[[1]])
    return(kp)
  }

  # Sort each group
  kp[] = lapply(kp, sortGroup)

  # Order groups by first element
  sortby = vapply(kp, function(g) g[1], 1L)
  if(anyDuplicated.default(sortby)) {
    sum = vapply(kp, function(g) sum(g), 1L)
    sortby = 1000L * sortby + sum
  }

  kp[.myorder(sortby, decreasing = TRUE)]
}


isDetailed = function(kp) {
  det = attr(kp, "detailed")
  !is.null(det) && det
}


kinpatReduce = function(kp, detailed = isDetailed(kp)) { # Remove empty groups

  # Remove empty blocks
  empty = lengths(kp) == 0
  if(any(empty))
    kp[empty] = NULL

  # If non-detailed, nothing more to do
  if(!detailed)
    return(kp)

  vec = unlist(kp, use.names = FALSE)

  # If no parental info, convert to non-detailed and return
  random = vec %% 10 == 0
  if(all(random)) {
    kp[] = lapply(kp, function(g) g %/% 10L)
    attr(kp, "detailed") = FALSE
    return(kp)
  }

  # Remove deterministic repeats in each block
  if(anyDuplicated.default(vec[!random])) {
    for(i in seq_along(kp)) {
      g = kp[[i]]
      dups = duplicated.default(g) & g %% 10 > 0
      if(any(dups))
        kp[[i]] = g[!dups]
    }
  }

  kp
}

# Replace first element with ancestor
kinpatRepl1 = function(x, anc) {
  x[[1]] = c(anc, x[[1]][-1])
  x
}

kinpatReplace = function(kp, id, rep1, rep2 = NULL) {
  g1 = kp[[1]]
  keep1 = g1 != id
  kp[[1]] = c(rep1, g1[keep1])

  if(!is.null(rep2)) {
    g2 = kp[[2]]
    keep2 = g2 != id
    kp[[2]] = c(rep2, g2[keep2])
  }
  kp
}


kinpatReplaceDet = function(kp, id, rep1, rep2 = NULL) {
  g1 = kp[[1]]
  keep1 = g1 %/% 10 != id
  kp[[1]] = c(rep1 * 10L, g1[keep1])

  if(!is.null(rep2)) {
    g2 = kp[[2]]
    keep2 = g2 %/% 10 != id
    kp[[2]] = c(rep2 * 10L, g2[keep2])
  }
  kp
}

#' @export
`[.kinpat` = function(x, idx, ...) {
  y = unclass(x)[idx]
  attributes(y) = attributes(x)
  y
}
