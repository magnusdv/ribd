context("Identity")

library(pedtools)
ident = function(x, ids, ...) condensedIdentity(x, ids, verbose=F, ...)

test_that("inbreeding coefficients of zero have no effect", {
  ans1 = c(2, 1, 4, 1, 4, 1, 7, 10, 2)/32
  ans2 = c(3, 2, 6, 1, 6, 1, 7, 6, 0)/32

  x = fullSibMating(1)

  founderInbreeding(x) = c('1'=0, '2'=0)
  expect_equal(ident(x, ids = 5:6), ans1)

  founderInbreeding(x) = c('1'=1, '2'=0)
  expect_equal(ident(x, ids = 5:6), ans2)

  founderInbreeding(x) = c('1'=0, '2'=1)
  expect_equal(ident(x, ids = 5:6), ans2)
})

test_that("sparse arrays give same answer", {
  ans1 = c(2, 1, 4, 1, 4, 1, 7, 10, 2)/32
  ans2 = c(1, 1, 2, 0, 2, 0, 2, 0, 0)/8

  x = fullSibMating(1)
  expect_equal(ident(x, 5:6, sparse=1), ans1)

  founderInbreeding(x, 1:2) = 1
  expect_equal(ident(x, 5:6, sparse=1), ans2)
})

test_that("PO with inbred parent gives correct answer with and without the full pedigree", {
  # Jacquard coeffs of father/son, when father has f=1/4.
  ans = c(0, 0, 1/4, 0, 0, 0, 0, 3/4, 0)

  # Full pedigree, creating father as a son of full sibs
  x1 = nuclearPed(2, sex=1:2)
  x1 = addChildren(x1, father=3, mother=4, nch=1)
  x1 = addSon(x1, 5)
  expect_equal(ident(x1, c(5,7)), ans)

  # Specifying the father's inbreeding coeff directly
  x2 = nuclearPed(1)
  founderInbreeding(x2, 1) = 1/4
  expect_equal(ident(x2, c(1,3)), ans)
})

test_that("founders with higher IDX than ids doesn't give error", {
  x = addChildren(nuclearPed(1), 3)

  founderInbreeding(x) = c('4'=1)
  expect_equal(ident(x, ids = 1:2), c(0,0,0,0,0,0,0,0,1))
})
