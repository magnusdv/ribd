stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

# Quick version of combn(., 2) for matrix output
comb2 = function(n, vec = length(n) > 1){
  if(vec) {
    v = n
    n = length(v)
  }
  if (n < 2)
    return(matrix(nrow = 0L, ncol = 2L))
  if (n == 2) {
    if(!vec) v = c(1L, 2L)
    return(`dim<-`(v, c(1L, 2L)))
  }
  if (n == 3) {
    x = c(1L, 1L, 2L, 2L, 3L, 3L)
    if(vec) x = v[x]
    return(`dim<-`(x, c(3L, 2L)))
  }

  x = rep.int(seq_len(n - 1), (n - 1):1)
  o = c(0, cumsum((n-2):1))
  y = seq_along(x) + 1 - o[x]

  if(vec)
    cbind(v[x], v[y], deparse.level = 0)
  else
    cbind(x, y, deparse.level = 0)
}

# A safer version of base::sample
safe_sample <- function(x, ...) x[sample.int(length(x), ...)]


# Fast setdiff
.mysetdiff = function(x, y) unique.default(x[match(x, y, 0L) == 0L])

# Fast intersection. NB: assumes no duplicates!
.myintersect = function(x, y) y[match(x, y, 0L)]

# Fast sorting of short vectors
.mysort = function(x, by = x, decreasing = FALSE) {
  len = length(x)
  if(length(by) != len)
    stop2("`x` and `by` must have the same length")

  if(len == 1)
    return(x)

  ord = .myorder(by, decreasing = decreasing)
  x[ord]
}

.myorder = function(x, decreasing = FALSE) {
  len = length(x)
  if(len > 3)
    return(order(x, decreasing = decreasing, method = "shell"))

  ord = switch(len,
    1,
    if(x[1] <= x[2]) 1:2 else 2:1,
    { # length 3
      a = x[1]; b = x[2]; c = x[3]
      ab = a <= b; ac = a <= c; bc = b <= c

      if(ab) {
        if(bc) c(1L,2L,3L)
        else if(ac) c(1L,3L,2L)
        else c(3L,1L,2L)
      }
      else {
        if(ac) c(2L,1L,3L)
        else if(bc) c(2L,3L,1L)
        else c(3L,2L,1L)
      }
    }
  )

  if(decreasing)
    ord = rev.default(ord)

  ord
}

# List pairs of IDs
.idPairs = function(v, self = FALSE, as = "character") {
  v = switch(as, character = as.character(v), integer = as.integer(v))
  n = length(v)
  sq = seq_along(v)

  # 2-column matrix
  if(self)
    mat = cbind(v[rep.int(sq, times = n:1)],
                v[unlist(lapply(sq, function(i) i:n))])
  else
    mat = comb2(v, vec = TRUE)

  lapply(1:nrow(mat), function(i) mat[i, ])
}


# Fast version of expand.grid
fast.grid = function(argslist, as.list = FALSE) {
  nargs = length(argslist)
  orep = nr = prod(lengths(argslist))
  if (nargs == 0L || nr == 0L)
    return(if(as.list) list() else matrix(ncol = 0, nrow = 0))

  rep.fac = 1L
  res = NULL
  for (x in argslist) {
    nx = length(x)
    orep = orep/nx
    res = c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)])  #this is res[, i]
    rep.fac = rep.fac * nx
  }
  dim(res) = c(nr, nargs)
  if (as.list)
    res = lapply(seq_len(nr), function(r) res[r, ])
  res
}

# Convert jacquard matrix to data frame
jmat2df = function(jmat, pairs, labs = NULL) {
  ids = unlist(pairs, use.names = FALSE)
  if(!is.null(labs))
    ids = labs[ids]
  dim(ids) = c(2, length(pairs))

  n = ncol(jmat)
  res = c(list(ids[1,], ids[2,]), lapply(1:n, function(i) jmat[,i]))
  names(res) = c("id1", "id2", paste0(if(n == 9) "D" else "d", 1:n))

  # Inspired by "quickdf" by Hadley W
  class(res) = "data.frame"
  attr(res, "row.names") = .set_row_names(length(pairs))

  res
}


# Add parents to selected founders
addFounderParents = function(x, ids, Xchrom = FALSE) {
  ids = unique.default(as.character(ids))

  idsFou = .myintersect(ids, founders(x))
  if(!length(idsFou))
    return(x)

  # Founder inbreeding?
  fouInb = founderInbreeding(x, idsFou, chromType = if(Xchrom) "X" else "autosomal")
  if(any(fouInb > 0))
    stop2("Cannot add parents to inbred founders")

  # Add parents
  for(id in idsFou)
    x = addParents(x, id, verbose = FALSE)

  x
}
