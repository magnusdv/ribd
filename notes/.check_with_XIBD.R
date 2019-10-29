#################################################################################
### Functions for comparing RIBD with (non-exported) methods in the XIBD package.
#################################################################################


compare_with_XIBD = function(x, ids, j) {
  message("Comparison with `XIBD` package: ", appendLF = FALSE)

  if(hasInbredFounders(x)) {
    message("skipped. (Pedigree has inbred founders.)")
    return()
  }

  jj = xibd(x, ids)
  if(is.null(jj)) {
    message("Install `XIBD` or use `checkAnswer = FALSE` to avoid this message.")
    return()
  }

  if(identical(j, jj))
    message("OK!")
  else if(isTRUE(all.equal(j,jj)))
    message("all.equal() OK, but not identical()")
  else {
    message("*** MISMATCH! ***")
    message("IDS: ", toString(ids))
    rbind(`XIBD:` = jj, `ribd:` = j)
  }
}



#' @importFrom utils capture.output
xibd = function(x, ids) {
  if(!requireNamespace("XIBD", quietly = TRUE)){
    message("Package `XIBD` is not installed.")
    return()
  }

  sex1 = getSex(x, ids[1])
  sex2 = getSex(x, ids[2])
  ids_int = internalID(x, ids)

  pedm = as.matrix(x)[, 1:4]
  colnames(pedm) = c("iid","pid", "mid","sex")

  capture.output(
    deltas <- lapply(1:9, function(i) XIBD:::delta(pedm, i, ids_int[1], ids_int[2]))
  )

  # If any of ids are male, the `deltas` list contains NULL entries
  jx = rep(NA_real_, 9)

  if(sex1 == 1 && sex2 == 1)
    jx[1:2] = unlist(deltas)
  else if (sex1 == 1 && sex2 == 2)
    jx[1:4] = unlist(deltas)
  else if (sex1 == 2 && sex2 == 1)
    jx[c(1,2,5,6)] = unlist(deltas)
  else
    jx[] = unlist(deltas)

  jx
}
