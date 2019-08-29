context("Two-locus identity coefficients")

tlid = function(x, ids, rho = 0.25, ...) {
  twoLocusIdentity(x, ids, rho = rho, ...)
}

test_that("two-loc-identity agrees with theory in FS with selfing", {
  x = ped(1:3, fid = c(0,1,1), mid = c(0,1,1), sex = c(0,1,1))
  J = condensedIdentity(x, 2:3)

  J2 = tlid(x, 2:3, rho = 0.25)
  expect_equal(rowSums(J2), J, check.attributes = F)
  expect_equal(colSums(J2), J, check.attributes = F)

  # Theory
  FSs = function(r) {
    rb = 1-r
    M = matrix(0, nrow = 9, ncol = 9, dimnames = rep(list(paste0("D", 1:9)), 2))
    M[1,1] = 1/8*(r^4 + rb^4)
    M[2,1] = M[1,2] = 1/4*rb^2*r^2
    M[3,1] = M[1,3] = 1/4*rb*r*(r^2 + rb^2)
    M[7,1] = M[1,7] = 1/2*rb^2*r^2

    M[5,] = M[,5] = M[3,]
    M[6,] = M[,6] = M[4,]
    M
  }

  expect_identical(J2[1, ], FSs(0.25)[1, ]) # TODO: finish this
})


test_that("two-loc-identity agrees with theory in full sib mating", {
  x = fullSibMating(1)
  J = condensedIdentity(x, 5:6)

  J2 = tlid(x, 5:6, rho = 0.25)
  expect_equal(rowSums(J2), J, check.attributes = F)
  expect_equal(colSums(J2), J, check.attributes = F)

  # Theory
  # TODO
})
