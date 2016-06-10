RSCRIPT = Rscript
TEX = pdflatex
PRINT = @echo -e "\e[1;34mBuilding $@\e[0m"
SOURCE_CPP_ALL := $(wildcard src/*.cpp)
SOURCE_CPP := $(filter-out src/RcppExports.cpp, $(SOURCE_CPP_ALL))
SOURCE_R_ALL := $(wildcard R/*.R)
SOURCE_R := $(filter-out R/RcppExports.R, $(SOURCE_R_ALL))

.PHONY: all test doc check_install rcpp vignette

all: install

test: install
	$(PRINT)
	R --vanilla -e 'devtools::test()'

install: doc
	$(PRINT)
	R --vanilla -e 'devtools::install()'

doc: R/RcppExports.R
	$(PRINT)
	R --vanilla -e 'devtools::document(roclets=c("rd", "collate", "namespace"))'

rcpp: R/RcppExports.R src/RcppExports.cpp

R/RcppExports.R src/RcppExports.cpp: $(SOURCE_CPP)
	$(PRINT)
	@echo -e $(SOURCE_CPP)
	R --vanilla -e 'Rcpp::compileAttributes()'
	touch R/RcppExports.R # because compileAttributes have is own makefile system
	touch src/RcppExports.cpp
