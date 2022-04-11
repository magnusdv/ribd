


  # non-inbred, X, female
  expect_equal(gk3_X(sf, c(1,1,1)), 0.25)
  expect_equal(gk4_X(sf, c(1,1,1,1)), 0.125)
  expect_equal(gk22_X(sf, c(1,1,1,1)), 0.25)

  # inbred, autosomal
  founderInbreeding(sf, 1) = ff = .25
  expect_equal(gk3(sf, c(1,1,1)), ff*1 +(1-ff)*0.25)
  expect_equal(gk4(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.125)
  expect_equal(gk22(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.25)

  # inbred, X, female
  founderInbreeding(sf, 1, chromType = "x") = ff = .25
  expect_equal(gk3_X(sf, c(1,1,1)), ff*1 +(1-ff)*0.25)
  expect_equal(gk4_X(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.125)
  expect_equal(gk22_X(sf, c(1,1,1,1)), ff*1 +(1-ff)*0.25)
})

