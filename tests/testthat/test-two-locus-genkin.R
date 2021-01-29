
tlgk = function(x, loc1, loc2 = "", rho = 0.25, ...)
  twoLocusGeneralisedKinship(x, loc1, loc2, rho = rho, ...)

test_that("two-loc generalised kinship, trivial examples", {
  x = nuclearPed(1)

  # s + t + u = 0
  expect_equal(tlgk(x, "3"), 1)
  expect_equal(tlgk(x, "3 = 1"), 1/4)
  expect_equal(tlgk(x, "3, 1"), 3/4)

  expect_equal(tlgk(x, "3>1 = 3>2 = 3>3"), 1/4)
  expect_equal(tlgk(x, "1 = 3>1 = 3>2 = 3>3"), 1/16)
  expect_equal(tlgk(x, "1, 3>1 = 3>2 = 3>3"), 3/16)

  # r > 0, s > 0, t = u = 0
  expect_equal(tlgk(x, "3>1, 3>2"), 1/2)
  expect_equal(tlgk(x, loc1="", loc2="3>1, 3>2"), 1/2)
  expect_equal(tlgk(x, "3>1 = 3>2, 3>3"), 1/4)  # r=2, s=1
  expect_equal(tlgk(x, "3>1, 3>2 = 3>3"), 1/4)  # r=1, s=2
  expect_equal(tlgk(x, "", "3>1 = 3>2, 3>3"), 1/4) # reverse loci
  expect_equal(tlgk(x, "", "3>1, 3>2 = 3>3"), 1/4)  # reverse loci

  # r > 0, t > 0, s = u = 0
  expect_equal(tlgk(x, "3>1", "3>1"), 1)
  expect_equal(tlgk(x, "3>1 = 3>2", "3>1 = 3>2"), 5/16)
  expect_equal(tlgk(x, "3>1 = 3>2", "3>1 = 3>3"), 1/4)

  # r > 0, s > 0, t > 0, u = 0
  expect_equal(tlgk(x, "3>1, 3>2", "3>1 = 3>2"), 3/16)
  expect_equal(tlgk(x, "3>1 = 3>2, 3>3", "3>1 = 3>2 = 3>3"), 3/64)

  # r > 0, s > 0, t > 0, u > 0
  expect_equal(tlgk(x, "3>1, 3>2", "3>1, 3>2"), 5/16)
})

test_that("two-loc gen kinship in inbred family", {
  x = halfCousinPed(0, child = T)
  rho = 0.25

  # Exact
  R = rho^2 + (1-rho)^2
  f = 1/8
  f11 = 1/8 * (1-rho)^2 * R
  f01 = f - f11
  f00 = 1 + f11 - 2*f

  # two-locus self kinship computed in two different ways
  a = tlgk(x, "6>1 = 6>2", "6>1 = 6>2", rho = rho)
  b = twoLocusKinship(x, c(6, 6), rho = rho)
  expect_equal(a, b)

  # two-locus kinship phi_00 (case r,s,t,u > 0)
  phi00 = tlgk(x, "6>1, 6>2", "6>1, 6>2", rho = rho)
  expect_equal(phi00, f00 * .5 * R)

  # slightly different, with 2 groups in L1 and 1 in L2. (In other words: u = 0)
  d = tlgk(x, "6>1, 6>2", "6>1 = 6>2", rho = rho)
  expect_equal(d, f01 * .5 + f00 * rho * (1-rho))
})


test_that("two-loc generalised kinship, selfing pedigree", {
  x = selfingPed(1)
  r = 0.25
  expect_equal(tlgk(x, "1>1 = 2>1", "1>1 = 2>1"), .5*(1-r)*(r^2 + (1-r)^2) + r/4)
})
