context("Sampler")



test_that("Sample from gentive model", {

  data <- SampleGenoFromGenerativeModel(30,1000,3,3)

  data <- SampleGenoFromGenerativeModel(30,1000,3,3,
                                        Q.sampler = SampleDistDFromCenterQ(0.1),
                                        G.sampler = SampleDirichletG(),
                                        coord.sampler = SampleNormalClusterCoord())

  data <- SampleGenoFromGenerativeModel(30,1000,3,3,
                                        Q.sampler = SampleDistFromCenterDirichletQ(),
                                        G.sampler = SampleDirichletG(),
                                        coord.sampler = SampleNormalClusterCoord())

  data <- SampleGenoFromGenerativeModel(30,1000,3,3,
                                        Q.sampler = SampleGaussianProcessQ(CovFunctionSquaredExp(0.2)),
                                        G.sampler = SampleDirichletG(),
                                        coord.sampler = SampleNormalClusterCoord())


})
