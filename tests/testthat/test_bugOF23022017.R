library(tess3r)
library(testthat)
context("Bug found by O.F. the 23/02/2017")

test_that("order of color in plot", {

  skip("run and see plot")

  qk <- as.qmatrix(read.table("./data-test/bugOF23022017/q3.txt"))
  coordinates <- as.matrix(read.table("./data-test/bugOF23022017/coord.txt"))

  my.colors  <- CreatePalette(c("red","green","blue"),9)

  barplot(qk, border = NA, space = 0,
          xlab = "",
          ylab = "", col.palette = my.colors) -> bp

  axis(1, at = 1:nrow(qk), labels = coordinates[,1], las = 2, cex.axis = .3)

  ## Les couleurs assignées naivement avec which.max. Les rouges sont a G, les bleus a D.

  ## c'est le bon mapping:

  plot(coordinates, pch = 19, col = 1 + apply(qk, 1, which.max)); map(add=T)

  # Avec le plot, les rouges sont a D, les bleus a G.
  ## Plus de bug -> les couleurs dans le même ordre ;-D

  plot(qk, coordinates, method = "map.max",
       interpol = FieldsKrigModel(10),
       main = "Ancestry coefficients",
       xlab = "Longitude", ylab = "Latitude",
       resolution = c(300,300), cex = .4,
       col.palette = my.colors
  )


  plot(qk, coordinates, method = "map.all",
       interpol = FieldsKrigModel(10),
       main = "Ancestry coefficients",
       xlab = "Longitude", ylab = "Latitude",
       resolution = c(300,300), cex = .4,
       col.palette = my.colors
  )

  ggtess3Q(qk, coordinates,
           col.palette = my.colors)
})
