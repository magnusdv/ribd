#' @rdname gKinship
#' @export
gip = function(x, pattern, distinct = TRUE) {
  if(!is.ped(x))
    stop2("First argument must be a `ped` object")

  if(isGip(pattern))
    pattern = gip2list(pattern)

  if(!is.list(pattern))
    stop2("`pattern` must be a list")

  names(pattern) = NULL # to avoid problems when unlisting
  ids = unlist(pattern)
  idsInt = internalID(x, ids)

  # Deterministic?
  nms = names(ids)
  determ = !is.null(nms)
  if(determ) {
    nmsNum = match(nms, c("", "p", "m")) - 1L  # 0, 1 and 2
    if(anyNA(nmsNum))
      stop2("Names (to indicate deterministic sampling) must be 'p' or 'm': ",
            setdiff(nms, c("", "p", "m")))
  }

  # Group index
  grp = rep(sq <- seq_along(pattern), lengths(pattern))

  # Create gip object
  pat = lapply(sq, function(i) {
    idx = grp == i
    g = idsInt[idx]
    if(determ)
      g = 10L*g + nmsNum[idx]
    g
  })

  structure(pat, labs = labels(x), deterministic = determ, distinct = distinct, class = "gip")
}


gip2string = function(gp, deterministic = isDeterministic(gp)) {

  labs = attr(gp, "labs")

  grps = vapply(gp, function(g) {
    if(deterministic) {
      lb = labs[g %/% 10]
      par = g %% 10
      lb[par > 0] = paste(lb[par > 0], c("p", "m")[par], sep = ":")
    }
    else {
      lb = labs[g]
    }
    sprintf("(%s)", paste0(lb, collapse = ","))
  }, FUN.VALUE = "")

  blocksep = if(isDistinct(gp)) " / " else " & "
  paste0(grps, collapse = blocksep)
}

#' @export
print.gip = function(x, ...) {
  cat(gip2string(x), "\n", sep = "")
}

gip2list = function(gp) {
  labs = attr(gp, "labs")
  determ = isDeterministic(gp)

  lapply(gp, function(g) {
    if(determ) {
      s = labs[g %/% 10]
      names(s) = c("", "p", "m")[g %% 10 + 1]
      s
    }
    else
      labs[g]
  })
}

gipSort = function(gp) {

  # Function for sorting a single group
  sortGroup = function(g) .mysort(g, decreasing = TRUE)

  # Quick return if just one group
  if(length(gp) == 1) {
    gp[[1]] = sortGroup(gp[[1]])
    return(gp)
  }

  # Sort each group
  gp[] = lapply(gp, sortGroup)

  # Order groups by first element
  sortby = vapply(gp, function(g) g[1], 1L)
  if(anyDuplicated.default(sortby)) {
    sum = vapply(gp, function(g) sum(g), 1L)
    sortby = 1000L * sortby + sum
  }

  gp[.myorder(sortby, decreasing = TRUE)]
}


isDeterministic = function(gp) {
  det = attr(gp, "deterministic")
  !is.null(det) && det
}

isDistinct = function(gp) {
  dist = attr(gp, "distinct")
  !is.null(dist) && dist
}

isGip = function(x) {
  inherits(x, "gip")
}

# Remove empty groups
gipReduce = function(gp, deterministic = isDeterministic(gp)) {

  # Remove empty blocks
  empty = lengths(gp) == 0
  if(any(empty))
    gp[empty] = NULL

  # If non-determ, nothing more to do
  if(!deterministic)
    return(gp)

  vec = unlist(gp, use.names = FALSE)

  # If no parental info, convert to non-determ and return
  random = vec %% 10 == 0
  if(all(random)) {
    gp[] = lapply(gp, function(g) g %/% 10L)
    attr(gp, "deterministic") = FALSE
    return(gp)
  }

  # Remove deterministic repeats in each block
  if(anyDuplicated.default(vec[!random])) {
    for(i in seq_along(gp)) {
      g = gp[[i]]
      dups = duplicated.default(g) & g %% 10 > 0
      if(any(dups))
        gp[[i]] = g[!dups]
    }
  }

  gp
}

# Replace first element with ancestor
gipRepl1 = function(x, anc) {
  x[[1]] = c(anc, x[[1]][-1])
  x
}

gipReplace = function(gp, id, rep1, rep2 = NULL) {
  g1 = gp[[1]]
  keep1 = g1 != id
  gp[[1]] = c(rep1, g1[keep1])

  if(!is.null(rep2)) {
    g2 = gp[[2]]
    keep2 = g2 != id
    gp[[2]] = c(rep2, g2[keep2])
  }
  gp
}


gipReplaceDet = function(gp, id, rep1, rep2 = NULL) {
  g1 = gp[[1]]
  keep1 = g1 %/% 10 != id
  gp[[1]] = c(rep1 * 10L, g1[keep1])

  if(!is.null(rep2)) {
    g2 = gp[[2]]
    keep2 = g2 %/% 10 != id
    gp[[2]] = c(rep2 * 10L, g2[keep2])
  }
  gp
}

#' @export
`[.gip` = function(x, idx, ...) {
  y = unclass(x)[idx]
  attributes(y) = attributes(x)
  y
}
