library(TESS3enchoSen)
# create dataset for testing
data.for.test = sampleTESS2.3(50,
                     2000,
                     1,
                     3 ,
                     SampleDistDFromCenterQ(0.1))
devtools::use_data(data.for.test)

# A thaliana data set
at.geno = "~/PatatorHomeDir/Data/At/little_sample/Athaliana.geno"
at.coord = "~/PatatorHomeDir/Data/At/little_sample/Athaliana.coord"
data.at = list()
data.at$X = LEA::read.geno(at.geno)
data.at$coord = as.matrix(read.table(at.coord))
devtools::use_data(data.at)
