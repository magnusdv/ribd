context("Kinship")

test_that("ibd_kappa gives message and returns NA,NA,NA if inbreeding", {
  x = fullSibMating(1)
  expect_message(ibd_kappa(x, 5:6), "Inbred individuals: 5, 6")
  expect_identical(ibd_kappa(x, 5:6), rep(NA_real_, 3))
})

test_that("ibd_kappa allows for rounding errors in the detection of inbreeding", {
  x = nuclearPed(2)
  f1 = 0.1
  f2 = 0.2

  # Exact kappa for sibs
  f00 = (1-f1)*(1-f2)
  f10 = f1*(1-f2)
  f01 = (1-f1)*f2
  f11 = f1*f2
  kappa_exact = c(f00*0.25, f00*0.5+f10*0.5+f01*0.5, f00*0.25+f10*0.5+f01*0.5+f11)

  # Computed
  founderInbreeding(x, 1:2) = c(f1, f2)
  kappa_ribd = ibd_kappa(x, ids = 3:4)

  # Results are not identical (due to rounding), but should be *equal*
  expect_equal(kappa_ribd, kappa_exact)
})
