context("genotype data functions")


test_that("test ComputeXBin", {
  M = matrix(as.integer(c(1,0,2,0,1,0,1,0,1)),3,3)
  MBin <- ComputeXBin(M, 2)
  expect_equal( mean( ComputeXFromXBin(MBin,2) - M), 0)
})
