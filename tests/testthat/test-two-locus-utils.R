context("Two-locus utilities")

test_that("two-loc kinship class replacment works", {
  x = linearPed(2)
  kin = kin2L(x, "3>5 = 4>5, 3>-1", locus2 = "3>5 = 3>-1")

  kin2 = kinReplace(kin, id=3, loc1Rep = list(from1 = 1, to1 = 301, from2 = 2, to2 = 302))
  expect_identical(kin2, kin2L(x, "1>301 = 4>5, 2>302", locus2 = "3>5 = 3>-1"))

})
