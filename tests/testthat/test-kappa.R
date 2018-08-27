context("Kinship")

library(pedtools)

test_that("ibd_kappa gives message and returns NA,NA,NA if inbreeding", {
  x = fullSibMating(1)
  expect_message(ibd_kappa(x, 5:6), "Inbred individuals: 5, 6")
  expect_identical(ibd_kappa(x, 5:6), rep(NA_real_, 3))
})

