wl = function(x, ...) gKinship(x, list(...), method = "WL")
ls = function(x, ...) gKinship(x, list(...), method = "LS")
gc = function(x, ...) gKinship(x, list(...), distinct = F, method = "GC")

test_that("Random distinct generalised kinships (WL) in singleton", {
  s = singleton(1)
  expect_equal(wl(s, c(1)), 1)
  expect_equal(wl(s, c(1,1)), 1/2)
  expect_equal(wl(s, 1,1), 1/2)
  expect_equal(wl(s, c(1,1,1,1)), 1/8)
  expect_equal(wl(s, c(1,1), c(1,1)), 1/8)
  expect_equal(wl(s, 1,1,1), 0)
})

test_that("Mixed distinct generalised kinships (LS) in singleton", {
  s = singleton(1)
  expect_equal(ls(s, c(1)), 1)
  expect_equal(ls(s, c(p=1,1)), 1/2)
  expect_equal(ls(s, c(m=1,m=1,1)), 1/2)
  expect_equal(ls(s, c(p=1,m=1,1)), 0)
  expect_equal(ls(s, 1, c(p=1)), 1/2)
  expect_equal(ls(s, 1, c(p=1, 1)), 1/4)
  expect_equal(ls(s, 1, c(p=1, m=1)), 0)
  expect_equal(ls(s, c(p=1), c(p=1)), 0)
  expect_equal(ls(s, c(p=1), c(m=1)), 1)
  expect_equal(ls(s, c(p=1,1), c(m=1)), 1/2)
  expect_equal(ls(s, c(p=1), c(m=1), 1), 0)
  expect_equal(ls(s, c(p=1), 1, 1), 0)
})

test_that("Deterministic nondistinct generalised kinships (GC) in singleton", {
  s = singleton(1)
  expect_equal(gc(s, c(p=1)), 1)
  expect_equal(gc(s, c(p=1,p=1)), 1)
  expect_equal(gc(s, c(p=1,m=1)), 0)
  expect_equal(gc(s, c(p=1,p=1), c(m=1,m=1)), 1)
})

test_that("Random distinct generalised kinships (WL) in selfingPed", {
  s = selfingPed(1)
  expect_equal(wl(s, c(2)), 1)
  expect_equal(wl(s, c(2,2)), 3/4)
  expect_equal(wl(s, 2,2), 1/4)
  expect_equal(wl(s, c(2,2,2)), .5 + .5*.25)
  expect_equal(wl(s, c(2,2), c(2,2)), .5 * 1/8)
  expect_equal(wl(s, c(1,2)), .5)
  expect_equal(wl(s, 1,2), .5)
  expect_equal(wl(s, 2,2,2), 0)
})

test_that("Mixed distinct generalised kinships (LS) in selfingPed", {
  s = selfingPed(1)
  f = 0.5
  expect_equal(ls(s, c(p=2)), 1)
  expect_equal(ls(s, c(p=2,2)), .75)
  expect_equal(ls(s, c(p=2,m=2)), .5)
  expect_equal(ls(s, c(p=2,m=2,2)), .5)
  expect_equal(ls(s, c(p=2),2), .25)
  expect_equal(ls(s, c(p=2),2,2), 0)
  expect_equal(ls(s, c(p=2),c(m=2)), .5)
  expect_equal(ls(s, c(p=2),c(p=2)), 0)
  expect_equal(ls(s, c(p=2),c(2,2)), 1/8)
})

