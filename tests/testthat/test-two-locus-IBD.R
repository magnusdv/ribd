
tlibd = function(x, ids, rho = 0.25, ...) {
  twoLocusIBD(x, ids, rho = rho, ...)
}

test_that("two-loc-IBD of PO agree with theory", {
  x = nuclearPed(1)

  # Theory
  PO = matrix(0, nrow = 3, ncol = 3, dimnames = rep(list(c("ibd0", "ibd1", "ibd2")), 2))
  PO[2,2] = 1

  expect_identical(tlibd(x, c(1,3)), PO)
  expect_identical(tlibd(x, c(1,3), rho = 0), PO)
  expect_identical(tlibd(x, c(1,3), rho = 0.5), PO)

  # Theory - detailed
  POdet = function(rho)
    c(k11.cc = 1-rho, k11.ct = 0, k11.tc = rho, k11.tt = 0)

  expect_equal(tlibd(x, c(1,3), coef = "k11", detailed = T, rho = 0),    POdet(0))
  expect_equal(tlibd(x, c(1,3), coef = "k11", detailed = T, rho = 0.25), POdet(0.25))
  expect_equal(tlibd(x, c(1,3), coef = "k11", detailed = T, rho = 0.5),  POdet(0.5))
})


test_that("two-loc-IBD of full sibs agree with theory", {
  x = nuclearPed(2)

  # Theory
  FS = function(rho) {
    R = rho^2 + (1-rho)^2
    nms = c("ibd0", "ibd1", "ibd2")
    m = matrix(0, nrow = 3, ncol = 3, dimnames = list(nms, nms))
    m[1,1] = m[3,3] = 0.25 *R^2
    m[2,1] = m[1,2] = 0.5 * R * (1-R)
    m[3,1] = m[1,3] = 0.25 * (1-R)^2
    m[2,2] = 0.5 * (1 - 2 * R * (1-R))
    m[3,2] = m[2,3] = 0.5 * R * (1-R)
    m
  }

  expect_identical(tlibd(x, 3:4, rho = 0),    FS(rho = 0))
  expect_identical(tlibd(x, 3:4, rho = 0.25), FS(rho = 0.25))
  expect_identical(tlibd(x, 3:4, rho = 0.5),  FS(rho = 0.5))

  # Theory - detailed
  FSdet_11 = function(r) {
    rb = 1-r
    R = r^2 + rb^2
    c(k11.cc = 0.5*R^2, k11.ct = 0, k11.tc = 0, k11.tt = 2*rb^2*r^2)
  }

  expect_identical(tlibd(x, 3:4, rho = 0, coef = "k11", detailed = T), FSdet_11(0))
  expect_identical(tlibd(x, 3:4, rho = 0.25, coef = "k11", detailed = T), FSdet_11(0.25))
  expect_identical(tlibd(x, 3:4, rho = 0.5, coef = "k11", detailed = T), FSdet_11(0.5))
})


test_that("two-loc-IBD of G with selfing, agrees with theory", {
  x = addSon(selfingPed(1), 2)

  # Theory -detailed
  Gs = function(rho) {
    rb = 1-rho
    c(k11.cc = rb^2 + rho/2, k11.ct = 0, k11.tc = rho*rb + rho/2, k11.tt = 0)
  }

  expect_identical(tlibd(x, c(1,4), rho = 0, coef = "k11"),    c(k11 = 1))
  expect_identical(tlibd(x, c(1,4), rho = 0.25, coef = "k11"), c(k11 = 1))
  expect_identical(tlibd(x, c(1,4), rho = 0.5, coef = "k11"),  c(k11 = 1))

  expect_identical(tlibd(x, c(1,4), rho = 0, coef = "k11", detailed = T),    Gs(0))
  expect_identical(tlibd(x, c(1,4), rho = 0.25, coef = "k11", detailed = T), Gs(0.25))
  expect_identical(tlibd(x, c(1,4), rho = 0.5, coef = "k11", detailed = T),  Gs(0.5))

})
