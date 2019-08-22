context("Two-locus kinship")

test_that("two-loc kinship is correct in selfing ped", {
  x = selfingPed(1)
  r = c(0, 0.25, 0.5)
  expect_identical(twoLocusKinship(x, ids = 1:2, rho = r),
                   .5*(1-r)*(r^2 + (1-r)^2) + .25*r)
})

test_that("two-loc kinship is correct in inbred example", {
  x = addChildren(nuclearPed(1, 2), 1, 3)
  phi11 = twoLocusKinship(x, ids = c(1,4), rho = 0.25)
  expect_identical(round(phi11, 5), 0.19238)
})
