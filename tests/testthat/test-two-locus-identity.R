
tlid = function(x, ids, rho = 0.25, ...) {
  twoLocusIdentity(x, ids, rho = rho, ...)
}

test_that("two-loc-identity agrees with theory in simple selfing ped", {
  x = selfingPed(1)
  J = condensedIdentity(x, 1:2)

  J2 = tlid(x, 1:2, rho = 0.25)
  expect_equal(rowSums(J2), J, check.attributes = F)
  expect_equal(colSums(J2), J, check.attributes = F)

  # Theory
  POs = function(r) {
    rb = 1-r
    M = matrix(0, nrow = 9, ncol = 9, dimnames = rep(list(paste0("D", 1:9)), 2))
    M[5,5] = M[7,7] = 1/2*(r^2 + rb^2)
    M[5,7] = M[7,5] = rb*r
    M
  }

  expect_identical(J2, POs(0.25))
})


test_that("two-loc-identity agrees with theory in FS with selfing", {
  x = ped(1:3, fid = c(0,1,1), mid = c(0,1,1), sex = c(0,1,1))
  J = condensedIdentity(x, 2:3)

  J2 = tlid(x, 2:3, rho = 0.25)
  expect_equal(rowSums(J2), J, check.attributes = F)
  expect_equal(colSums(J2), J, check.attributes = F)

  # Theory
  FSs = function(r) {
    rb = 1-r
    R = r^2 + rb^2
    M = matrix(0, nrow = 9, ncol = 9, dimnames = rep(list(paste0("D", 1:9)), 2))
    M[1,1] = 1/8*(r^4 + rb^4)
    M[1,2] = M[2,1] = 1/4*rb^2*r^2
    M[1,3] = M[3,1] = 1/4*rb*r*R
    M[1,5] = M[5,1] = 1/4*rb*r*R
    M[1,7] = M[7,1] = 1/2*rb^2*r^2

    M[2,2] = 1/8*(rb^4 + r^4)
    M[2,3] = M[3,2] = 1/4*rb*r*R
    M[2,5] = M[5,2] = 1/4*rb*r*R
    M[2,7] = M[7,2] = 1/2*rb^2*r^2

    M[3,3] = 1/4*R^2
    M[3,5] = M[5,3] = rb^2*r^2
    M[3,7] = M[7,3] = 1/2*rb*r*R

    M[5,5] = 1/4*R^2
    M[5,7] = M[7,5] = 1/2*rb*r*R

    M[7,7] = 1/4*R^2

    M[4,] = M[,4] = M[6,] = M[,6] = M[8,] = M[,8] = M[9,] = M[,9] = 0
    M
  }

  expect_equal(J2, FSs(0.25))
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
