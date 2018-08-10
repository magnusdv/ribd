context("Kinship")

library(pedtools)

test_that("kinship coefficients are the same with ribd and kinship2", {
  x = randomPed(4, founders=2)
  expect_identical(ibd_kinship(x), kinship2_kinship(x))
})


test_that("kinship coefficients with inbred founders are correct", {
  x = nuclearPed(fa="fa", mo="mo", child="boy")
  x = reorderPed(x, 3:1) # make it harder

  # make father a child of full sibs
  y1 = relabel(cousinsPed(0, child=T), c(101:104, "fa"))
  x_big = mergePed(x, y1)

  # make mother a child of half sibs
  y2 = relabel(halfCousinsPed(0, child=T), c(201:205, "mo"))
  y2 = swapSex(y2, "mo")
  x_big = mergePed(x_big, y2)

  ans1 = ibd_kinship(x_big)[x$LABELS, x$LABELS]

  # With founder inbreeding:
  founder_inbreeding(x) = c(fa=1/4, mo=1/8)
  ans2 = ibd_kinship(x)

  expect_identical(ans1, ans2)
})
