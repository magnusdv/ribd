context("Constructing pedigrees from kappa")

testConstruction = function(k0, k2) {
  k = c(k0, 1-k0-k2, k2)
  x = constructPedigree(k, verbose = FALSE)
  k2 = kappaIBD(x, leaves(x))
  all.equal(k, k2)
}

test_that("constructPedigree() reproduces kappa", {
  expect_true(testConstruction(0, 0))
  expect_true(testConstruction(0, 1))
  expect_true(testConstruction(0, 0.5))
  expect_true(testConstruction(0.5, 0))
  expect_true(testConstruction(0.01, 0.01))
  expect_true(testConstruction(0.25, 0.25))
  expect_true(testConstruction(17/32, 1/32))
})
