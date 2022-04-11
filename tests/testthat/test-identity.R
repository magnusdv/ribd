
ident = function(x, ids, method = c("K", "WL", "LS", "GC"), skip = character(), ...) {
  j0 = condensedIdentity(x, ids, verbose = F, ...)
  for(m in setdiff(method, skip)) {
    j = identityCoefs(x, ids, method = m, ...)
    # print(m); print(j)
    if(!isTRUE(all.equal(j, j0)))
      stop2(sprintf("Wrong answer with method %s\nj0: %s\n%s: %s", m, toString(j0), m, toString(j)))
  }
  j0
}

test_that("jacquard coeffs are correct in fullsibmating", {
  x = fullSibMating(1)

  ans12 = c(0,0,0,0,0,0,0,0,1)
  ans15 = c(0,0,0,0,1,1,1,4,1)/8
  ans35 = c(0,0,0,0,1,0,1,2,0)/4
  ans56 = c(2, 1, 4, 1, 4, 1, 7, 10, 2)/32

  # Normal ordering
  expect_equal(ident(x, ids = 1:2), ans12)
  expect_equal(ident(x, ids = c(1,5)), ans15)
  expect_equal(ident(x, ids = c(3,5)), ans35)
  expect_equal(ident(x, ids = 5:6), ans56)

  # Reordered
  y = reorderPed(x, 6:1)
  expect_equal(ident(y, ids = 1:2), ans12)
  expect_equal(ident(y, ids = c(1,5)), ans15)
  expect_equal(ident(y, ids = c(3,5)), ans35)
  expect_equal(ident(y, ids = 5:6), ans56)
})


test_that("jacquard with founder inbreeding, fullsibmating", {
  ans1 = c(3, 2, 6, 1, 6, 1, 7, 6, 0)/32
  ans2 = c(1,1,2,0,2,0,2,0,0)/8

  x = fullSibMating(1)

  founderInbreeding(x) = c('1'=1, '2'=0)
  expect_equal(ident(x, ids = 5:6), ans1)

  founderInbreeding(x) = c('1'=1, '2'=1)
  expect_equal(ident(x, ids = 5:6), ans2)
})

test_that("sparse arrays give same answer", {
  ans1 = c(2, 1, 4, 1, 4, 1, 7, 10, 2)/32
  ans2 = c(1, 1, 2, 0, 2, 0, 2, 0, 0)/8

  x = fullSibMating(1)
  expect_equal(ident(x, 5:6, method = "K", sparse=1), ans1)

  founderInbreeding(x, 1:2) = 1
  expect_equal(ident(x, 5:6, method = "K", sparse=1), ans2)
})

test_that("PO with inbred parent, with and without full pedigree", {
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
  expect_equal(ident(x2, c(1,3), skip = "LS"), ans)
})



# X chromosomal -----------------------------------------------------------

identX = function(x, ids, method = c("K", "GC"), skip = character(), ...) {
  j0 = condensedIdentityX(x, ids, verbose = F, ...)
  for(m in setdiff(method, skip)) {
    j = identityCoefs(x, ids, Xchrom = T, method = m, ...)
    # print(m); print(j)
    if(!isTRUE(all.equal(j, j0)))
      stop2(sprintf("Wrong answer with method %s\nj0: %s\n%s: %s", m, toString(j0), m, toString(j)))
  }
  j0
}


test_that("jacquard-X in fullsibmating", {
  x = fullSibMating(1)

  ans12 = c(0,0,0,1,rep(NA,5))
  ans15 = c(0.5,0.5,rep(NA, 7))
  ans16 = c(0,1,2,1,rep(NA,5))/4
  ans26 = c(0,0,0,0,1,0,1,2,0)/4
  ans35 = c(1,3,rep(NA,7))/4
  ans45 = c(0,0,NA,NA,1,0,NA,NA,NA)
  ans46 = c(0,0,0,0,1,0,1,2,0)/4
  ans56 = c(1,1,4,2,rep(NA,5))/8

  # Normal ordering
  expect_equal(identX(x, ids = 1:2), ans12)
  expect_equal(identX(x, ids = c(1,5)), ans15)
  expect_equal(identX(x, ids = c(1,6)), ans16)
  expect_equal(identX(x, ids = c(2,6)), ans26)
  expect_equal(identX(x, ids = c(3,5)), ans35)
  expect_equal(identX(x, ids = c(4,5)), ans45)
  expect_equal(identX(x, ids = c(4,6)), ans46)
  expect_equal(identX(x, ids = 5:6), ans56)
})


test_that("jacquard-X with X-founder-inbreeding, fullsibmating", {
  x = fullSibMating(1)
  founderInbreeding(x, chromType = "x") = c('2'=1)

  expect_equal(identX(x, ids = c(2,6)), c(.5,0,.5,rep(0,6)))
  expect_equal(identX(x, ids = 5:6), c(1,1,2,0,rep(NA,5))/4)
})
