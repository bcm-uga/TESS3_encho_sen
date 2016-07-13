FROM rocker/hadleyverse

MAINTAINER cayek "https://github.com/cayek"

################################################################################
# install dependencies
RUN R -e 'install.packages(c("raster", \
"maps", \
"permute"))'

################################################################################
# install tess3
RUN R -e 'devtools::install_github("BioShock38/TESS3_encho_sen@master")'
