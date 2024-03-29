
test_that("kinship coefficients are the same with ribd and kinship2", {
  x = randomPed(10, seed=1234)
  expect_identical(kinship(x), kinship2_kinship(x))
  expect_identical(kinship(x, Xchrom = T), kinship2_kinship(x, Xchrom = T))
})

test_that("kinship() gives same result with and without `ids`", {
  x = quadHalfFirstCousins()
  expect_identical(kinship(x, ids = 9:10), kinship(x)[9,10])
  expect_identical(kinship(x, ids = c(1,1)), kinship(x)[1,1])

  set.seed(123)
  y = reorderPed(x, sample(10, ))
  expect_identical(kinship(y, ids = 9:10), kinship(y)["9","10"])
  expect_identical(kinship(y, ids = c(1,1)), kinship(y)["1","1"])
})

test_that("kinship coefficients with inbred founders are correct", {
  x = nuclearPed(fa="fa", mo="mo", child="boy")
  x = reorderPed(x, 3:1) # make it harder

  # make father a child of full sibs
  y1 = relabel(cousinPed(0, child=T), c(101:104, "fa"))
  x_big = mergePed(x, y1)

  # make mother a child of half sibs
  y2 = relabel(halfCousinPed(0, child=T), c(201:205, "mo"))
  y2 = swapSex(y2, "mo")
  x_big = mergePed(x_big, y2)

  labs = labels(x)
  ans1 = kinship(x_big)[labs, labs]

  # With founder inbreeding:
  founderInbreeding(x) = c(fa=1/4, mo=1/8)
  ans2 = kinship(x)

  expect_identical(ans1, ans2)
})


test_that("inbred founders are detected in `inbreeding()`", {
  # No founder inbreeding
  x = nuclearPed(1)
  expect_identical(inbreeding(x), structure(c(0,0,0), names=1:3))

  # With founder inbreeding
  y = x
  founderInbreeding(y, 2) = 1
  expect_identical(inbreeding(y), structure(c(0,1,0), names=1:3))
})

test_that("inbreeding coefficients are correctly computed", {
  x = y = cousinPed(0, child = T)
  expect_identical(inbreeding(x), structure(c(0,0,0,0,1/4), names=1:5))

  # With founder inbreeding
  founderInbreeding(y, 1) = 1
  expect_identical(inbreeding(y), structure(c(1,0,0,0,3/8), names=1:5))
})

test_that("inbreeding() works in selfing pedigree", {
  s = selfingPed(1)
  expect_equal(inbreeding(s, ids = 1), 0)
  expect_equal(inbreeding(s, ids = 2), 0.5)
})

test_that("X-chrom inbreeding is computed correctly", {
  x = halfCousinPed(0, child = T)
  xMat = swapSex(x, 1)
  child = leaves(x)
  fou = commonAncestors(x, parents(x, child)) # robust

  # Always 1 for males
  expect_equal(inbreeding(x, ids = child, Xchrom = T), 1)
  expect_equal(inbreeding(xMat, ids = child, Xchrom = T), 1)

  # Female child
  x = swapSex(x, child)
  xMat = swapSex(xMat, child)
  expect_equal(inbreeding(x, ids = child, Xchrom = T), 0)
  expect_equal(inbreeding(xMat, ids = child, Xchrom = T), 0.25)

  # With founder inbreeding
  founderInbreeding(xMat, fou, chrom = "x") = 1
  expect_equal(inbreeding(xMat, id = 6, Xchrom = T), 0.5)
})


test_that("kinship() works in pedlist", {
  x1 = nuclearPed() |> setFounderInbreeding(1, value = 1) |> reorderPed(3:1)
  x2 = singleton(4)
  x3 = singleton("NN") |> setFounderInbreeding(value = 0.5)
  x = list(x1, x2, x3)

  expect_equal(kinship(x)[1:3, 1:3], kinship(x1))
  expect_equal(kinship(x, x1$ID), kinship(x1))
  expect_equal(kinship(x, 1:3), kinship(x1, ids = 1:3))
  expect_equal(kinship(x, 4), kinship(x2))
  expect_equal(kinship(x)[4, 4, drop = F], kinship(x2))

  ids = c("NN",3,2)
  expect_equal(kinship(x, ids), kinship(x)[ids, ids])

  ids = c("NN",4)
  expect_equal(kinship(x, ids, simplify = F), kinship(x)[ids, ids])

  expect_equal(kinship(x, 2:3), 0.25)
  expect_equal(kinship(x, 3:4), 0)

  expect_error(kinship(x, 4:5), "Unknown ID label: 5")
})


test_that("kinship(Xchrom = T) works in pedlist", {
  x1 = nuclearPed() |> setFounderInbreeding(1, value = 1) |> reorderPed(3:1)
  x2 = singleton(4)
  x3 = singleton("NN") |> setFounderInbreeding(value = 0.5)
  x = list(x1, x2, x3)

  kinX = function(...) kinship(..., Xchrom = T)

  expect_equal(kinX(x)[1:3, 1:3], kinX(x1))
  expect_equal(kinX(x, x1$ID), kinX(x1))
  expect_equal(kinX(x, 1:3), kinX(x1, ids = 1:3))
  expect_equal(kinX(x, 4), kinX(x2))
  expect_equal(kinX(x)[4, 4, drop = F], kinX(x2))

  ids = c("NN",3,2)
  expect_equal(kinX(x, ids), kinX(x)[ids, ids])

  ids = c("NN",4)
  expect_equal(kinX(x, ids, simplify = F), kinX(x)[ids, ids])

  expect_equal(kinX(x, 2:3), 0.5)
  expect_equal(kinX(x, 3:4), 0)

  expect_error(kinX(x, 4:5), "Unknown ID label: 5")
})

test_that("kinship() works in pedlist with duplicated IDs", {
  x1 = nuclearPed(children = 3:4)
  x2 = nuclearPed(children = 6:5)
  x = list(x1, x2)

  expect_equal(kinship(x, 3:4), 0.25)
  expect_equal(kinship(x, 4:5), 0)
  expect_equal(kinship(x, x1$ID, simplify = F), kinship(x1))
  expect_error(kinship(x), "ID label is not unique: 1")
})
