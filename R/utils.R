stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Quick version of combn(., 2)
comb2 = function(n, vec = length(n) > 1){
  if(vec) {
    v = n
    n = length(v)
  }
  if (n < 2)
    return(matrix(nrow = 0, ncol = 2))

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

# TODO: Move to pedtools
foundersFirst = function(x) {
  fou = founders(x, internal = TRUE)

  # Check if all foundders are already first
  if(length(fou) == max(fou))
    return(x)

  nonfou = nonfounders(x, internal = TRUE)
  reorderPed(x, neworder = c(fou, nonfou))
}

