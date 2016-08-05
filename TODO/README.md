# Olivier
- vignette (there are in the dir vignettes/quarantine/)
- example of function must compile (barplot.tess3Q)

# Kevin
- make an exemple/test with the minimal example :
``` R
obj <- tess3(matrix, coord, K = 1:6)
plot(obj) #plot l'erreur de CV, choix de K
q.matrix <- qmatrix(obj, K = 3) #get Q
barplot(q.matrix, K = 3)
plot(q.matrix, coord, K=3) #map les coeffs par krig
p.values = pvalue(obj, K=3) #get les pvaleurs
```
- plot of project : plot(obj) #plot l'erreur de CV, choix de K
- getter to return Q, G etc cad :
  - qmatrix(obj , K, run ="best")
  - pvalue(obj , K, run ="best")
  - ...
- option names must be consistent with name in article
- check function : control class : "matrix" %in% class(Q) etc...
- fix warning and note when running make check
- option chisq = TRUE (use chisq to compute pvalue)

# Flora
- ...

# Someone
- vignette for different data format files
