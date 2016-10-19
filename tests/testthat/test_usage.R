context("Usage")


test_that("Minimal example", {
  data("data.for.test", package = "tess3r")
  obj <- tess3(X = data.for.test$X,
               coord = data.for.test$coord,
               K = 1:6,
               ploidy = 1)
  plot(obj) #plot l'erreur de CV, choix de K
  q.matrix <- qmatrix(obj, K = 3) #get Q
  expect_equal(dim(q.matrix), c(data.for.test$n, 3))
  barplot(q.matrix)
  plot(q.matrix, data.for.test$coord) #map les coeffs par krig
  # FJ: previous line does not work for me with backgound=T because it's in the sea, adding following:
  plot(q.matrix, data.for.test$coord, background = FALSE)
  p.values = pvalue(obj, K=3) #get les pvaleurs
  expect_equal(dim(p.values), c(data.for.test$L, 1))
  piechart(q.matrix, data.for.test$coord) #map les coeffs par krig


  })
