
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



identity_merlin = function(x, ids = labels(x), self = FALSE, detailed = FALSE, execPath = "merlin", verbose = FALSE, cleanup = TRUE) {

  # Check merlin availability
  if(!nzchar(pth <- Sys.which(execPath))) {
    stop2("MERLIN not found. Use the `execPath` argument to supply the path to the executable")
  }

  x = prepPed(x)
  pairs = prepIds(x, ids, self = self, int = TRUE)

  # Add empty marker (seems to be needed for MERLIN to work)
  x = setMarkers(x, NULL)
  x = addMarker(x)

  prefix = "__ribd2merlin__"
  fin = paste(prefix, c("ped", "map", "dat"), sep = ".")
  fout = paste(prefix, "s15", sep = ".")

  if(cleanup)
    on.exit({unlink(c(fin, fout))})

  ### Generate MERLIN input files
  # ped file
  pedmatr = cbind(1, as.matrix(x, include.attrs = FALSE))
  write(t.default(pedmatr), file = fin[1], ncolumns = ncol(pedmatr))

  ### map file
  write(c(CHROM = 1, MARKER = NA, MB = NA), file = fin[2], ncolumns = 3)

  ### dat file
  write(c("M", NA), file = fin[3], ncolumns = 2)

  # Run MERLIN --extended
  commandArgs = c(execPath,
                  sprintf("-p %s.ped -d %s.dat -m %s.map --prefix %s", prefix, prefix, prefix, prefix),
                  sprintf("--bits %d", 2*pedsize(x) - length(founders(x))),
                  "--extended")
  command = paste(commandArgs, collapse = " ")

  if (verbose)
    cat("\nExecuting the following command:\n", paste0(commandArgs, collapse = "\n "), "\n", sep = "")

  mout = suppressWarnings(system(command, intern = TRUE))

  # Catch errors
  err = NULL
  if (any(fatal <- substr(mout, 1, 11) == "FATAL ERROR"))
    err = mout[which(fatal)[1]:length(mout)]
  else if (any(warn <- substr(mout, 2, 8) == "WARNING"))
    err = mout[which(warn)[1] + 0:5]

  if(!is.null(err))
    warning(paste0(err, collapse = "\n"), call. = FALSE)
  else if (verbose)
    cat("\nMERLIN run completed\n")

  # Load results file
  j = read.table(fout, header = TRUE)

  # Extracted wanted pairs/cols
  j[j$ID1 > j$ID2, 2:3] = j[j$ID1 > j$ID2, 3:2] # swap ids
  rownames(j) = paste(j$ID1, j$ID2, sep = "-")
  jsub = j[sapply(pairs, paste, collapse="-"), , drop = FALSE]

  # Prepare output
  labs = labels(x)
  res = data.frame(id1 = labs[jsub$ID1], id2 = labs[jsub$ID2], jsub[,-(1:4)])
  names(res)[-(1:2)] = paste0("d", seq_len(ncol(res) - 2)) # not sure if always 15
  rownames(res) = NULL

  if(!detailed)
    res = detailed2condensed(res)

  res
}


