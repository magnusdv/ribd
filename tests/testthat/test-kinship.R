context("Kinship")

test_that("kinship coefficients are the same with ribd and kinship2", {
  x = randomPed(4, founders=2)
  expect_identical(ibd_kinship(x), kinship2_kinship(x))
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
  ans1 = ibd_kinship(x_big)[labs, labs]

  # With founder inbreeding:
  founderInbreeding(x) = c(fa=1/4, mo=1/8)
  ans2 = ibd_kinship(x)

  expect_identical(ans1, ans2)
})


test_that("inbred founders are detected in `ibd_inbreeding()`", {
  # No founder inbreeding
  x = nuclearPed(1)
  expect_identical(ibd_inbreeding(x), structure(c(0,0,0), names=1:3))

  # With founder inbreeding
  y = x
  founderInbreeding(y, 2) = 1
  expect_identical(ibd_inbreeding(y), structure(c(0,1,0), names=1:3))
})

test_that("inbreeding coefficients are correctly computed", {
  x = y = cousinPed(0, child = T)
  expect_identical(ibd_inbreeding(x), structure(c(0,0,0,0,1/4), names=1:5))

  # With founder inbreeding
  founderInbreeding(y, 1) = 1
  expect_identical(ibd_inbreeding(y), structure(c(1,0,0,0,3/8), names=1:5))
})
