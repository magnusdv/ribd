context("Generalised kinships")

gk3 = function(x, ids, ...) generalisedKinship3(x, ids, verbose=F, ...)
gk4 = function(x, ids, ...) generalisedKinship4(x, ids, verbose=F, ...)
gk22 = function(x, ids, ...) generalisedKinship22(x, ids, verbose=F, ...)

gk3_X = function(x, ids, ...) generalisedKinship3(x, ids, chromType = "x", verbose=F, ...)
gk4_X = function(x, ids, ...) generalisedKinship4(x, ids, chromType = "x", verbose=F, ...)
gk22_X = function(x, ids, ...) generalisedKinship22(x, ids, chromType = "x", verbose=F, ...)

test_that("generalised kinship coefs of a singleton are correct", {
  sm = singleton(1, sex=1)
  sf = singleton(1, sex=2)

  # non-inbred, autosomal
  expect_equal(gk3(sm, c(1,1,1)), 0.25)
  expect_equal(gk4(sm, c(1,1,1,1)), 0.125)
  expect_equal(gk22(sm, c(1,1,1,1)), 0.25)

  # non-inbred, X, male
  expect_equal(gk3_X(sm, c(1,1,1)), 1)
  expect_equal(gk4_X(sm, c(1,1,1,1)), 1)
  expect_equal(gk22_X(sm, c(1,1,1,1)), 1)

  # non-inbred, X, female
  expect_equal(gk3_X(sf, c(1,1,1)), 0.25)
  expect_equal(gk4_X(sf, c(1,1,1,1)), 0.125)
  expect_equal(gk22_X(sf, c(1,1,1,1)), 0.25)

  # inbred, autosomal
  founderInbreeding(sf, 1) = ff = .25
  expect_equal(gk3(sf, c(1,1,1)), ff*1 +(1-ff)*0.25)
  expect_equal(gk4(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.125)
  expect_equal(gk22(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.25)

  # inbred, X, female
  founderInbreeding(sf, 1, chromType = "x") = ff = .25
  expect_equal(gk3_X(sf, c(1,1,1)), ff*1 +(1-ff)*0.25)
  expect_equal(gk4_X(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.125)
  expect_equal(gk22_X(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.25)
})

