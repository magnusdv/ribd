context("Two-locus utilities")

test_that("two-loc kinship class replacment works", {
  x = linearPed(2)
  kin = kin2L(x, "3>5 = 4>5, 3>-1", locus2 = "3>5")

  kin2 = kinRepl(kin, id=3, par1=1:2, gr1=1:2, par2=1, gr2 = 1)
  expect_identical(kin2, kin2L(x, "1>3 = 4>5, 2>3", locus2 = "1>3"))

})
