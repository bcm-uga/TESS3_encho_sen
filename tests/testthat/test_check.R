context("check")


test_that("cheks", {
  # X
  data("data.for.test", package = "tess3r")
  CheckX(data.for.test$X, 1)

  # W
  W <- ComputeHeatKernelWeight(coord = data.for.test$coord,
                               sigma = 2.0)
  CheckW(W)

  # W X
  CheckXW(data.for.test$X, W)

  # raise error
  # W
  W[1,2] <- 1000
  expect_error(CheckW(W), "W must be a squared symmetric matrix")
  # W X
  W <- matrix(0,1,1)
  expect_error(CheckXW(data.for.test$X, W), "*W must be of size*")
  # W X Coord
  W <- matrix(0, nrow(data.for.test$X), nrow(data.for.test$X))
  coord <- matrix(1,10,2)
  expect_error(CheckWCoord(W, coord), "Number of row of coord and W must be the same")

})
