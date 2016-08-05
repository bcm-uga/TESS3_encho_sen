FROM rocker/hadleyverse

MAINTAINER cayek "https://github.com/cayek"

################################################################################
# install dependencies
RUN install2.r --error \
    -r "http://cran.rstudio.com" \
    -r "http://www.bioconductor.org/packages/release/bioc" \
    testthat \
    sp \
    raster \
    knitr \
    rmarkdown \
    maps \
    permute \
    RColorBrewer \
    gtools \
    MASS \
    LEA \
    mapplots \
    gstat \
    sp \
    automap \
    fields \
    Rcpp \
    RcppEigen \
    quadprog

################################################################################
# install tess3
RUN R -e 'devtools::install_github("BioShock38/TESS3_encho_sen@master")'
