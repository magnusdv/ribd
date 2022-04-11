
identity_identity = function(x, ids, self = FALSE, verbose = FALSE) {
  if (!requireNamespace("identity", quietly = TRUE))
    stop2("Package `identity` must be installed for this function to work")

  x = prepPed(x)
  ped = cbind(1:pedsize(x), x$FIDX, x$MIDX)

  # This doesn't work due to bug in identity::identity.coefs():
  # pairs = prepIds(x, ids, self = self, int = TRUE)
  # pmat = do.call(rbind, pairs)
  # j = identity::identity.coefs(pmat, ped)

  idsInt = .mysort(internalID(x, unique.default(ids)))
  j = identity::identity.coefs(idsInt, ped)

  labs = labels(x)
  res = data.frame(id1 = labs[j[,1]], id2 = labs[j[,2]], j[, 3:11])
  names(res)[3:11] = paste0("D", 1:9)

  if(!self && nrow(res) > 1) {
    res = res[res$id1 != res$id2, , drop = FALSE]
    rownames(res) = NULL
  }

  res
}


#' @importFrom utils read.table write.table
identity_idcoefs = function(x, ids, self = FALSE, execPath = "idcoefs", verbose = FALSE, ram = 100, cleanup = TRUE) {

  # Check availability of `idcoefs` program
  if(!nzchar(pth <- Sys.which(execPath))) {
    stop2("Executable not found. Use the `execPath` argument to supply the path to the idcoefs executable")
  }

  x = prepPed(x)
  ped = cbind(1:pedsize(x), x$FIDX, x$MIDX)

  pairs = prepIds(x, ids, self = self, int = TRUE)
  pmat = do.call(rbind, pairs)

  # Clean now and possibly on exit
  clean = function() unlink(dir(pattern = "__ribd2idcoefs__"))
  clean()
  if(cleanup)
    on.exit(clean())

  # Write files
  write.table(ped, file = "__ribd2idcoefs__.ped",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(pmat, file = "__ribd2idcoefs__.sample",
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Run idCoefs
  command = paste(execPath,
                  "-p __ribd2idcoefs__.ped",
                  "-s __ribd2idcoefs__.sample",
                  "-o __ribd2idcoefs__.output",
                  "-r", ram)

  tryCatch({
    run = system(command, intern = TRUE)
    if(verbose)
      cat(run, sep = "\n")
    j = read.table("__ribd2idcoefs__.output", as.is = TRUE)
  },
    error = function(e) stop2(conditionMessage(e)),
  warning = function(e) stop2(conditionMessage(e))
  )

  # Parse results
  labs = labels(x)
  res = data.frame(id1 = labs[j[,1]], id2 = labs[j[,2]], j[, 3:11])
  names(res)[3:11] = paste0("D", 1:9)

  res
}



