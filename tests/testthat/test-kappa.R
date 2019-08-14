context("Kinship")

test_that("kappaIBD() deals sensibly with inbred individuals", {
  x = nuclearPed(1)
  founderInbreeding(x, 1:2) = 1
  # action = 0
  expect_identical(kappaIBD(x, 1:2, inbredAction = 0), rep(NA_real_, 3))

  # action = 1
  expect_identical(kappaIBD(x, 1:2, inbredAction = 1), rep(NA_real_, 3))
  expect_message(kappaIBD(x, 1:2, inbredAction = 1), "Warning:")

  # action = 2
  expect_error(kappaIBD(x, 1:2, inbredAction = 2), "Kappa coefficients are only")

})

test_that("kappa allows for rounding errors in the detection of inbreeding", {
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
  kappa_ribd = kappaIBD(x, ids = 3:4)

  # Results are not identical (due to rounding), but should be *equal*
  expect_equal(kappa_ribd, kappa_exact)
})
