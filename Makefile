RED = \033[0;31m
GREEN =  \033[0;32m
NC = \033[0m # No Color
PRINT = @echo "$(GREEN)Building $@ $(NC)"

################################################################################
# compile

SOURCE_CPP_ALL := $(wildcard src/*.cpp)
SOURCE_CPP := $(filter-out src/RcppExports.cpp, $(SOURCE_CPP_ALL))
SOURCE_R_ALL := $(wildcard R/*.R)
SOURCE_R := $(filter-out R/RcppExports.R, $(SOURCE_R_ALL))

.PHONY: all test doc check_install rcpp vignette check

all: install

test:
	$(PRINT)
	R --vanilla -e 'devtools::test()'

check:
	$(PRINT)
	R --vanilla -e 'devtools::check()'

install:
	$(PRINT)
	R --vanilla -e 'devtools::install()'

doc:
	$(PRINT)
	R --vanilla -e 'devtools::document(roclets=c("rd", "collate", "namespace"))'

rcpp: R/RcppExports.R src/RcppExports.cpp

R/RcppExports.R src/RcppExports.cpp: $(SOURCE_CPP)
	$(PRINT)
	@echo -e $(SOURCE_CPP)
	R --vanilla -e 'Rcpp::compileAttributes()'
	touch R/RcppExports.R # because compileAttributes have is own makefile system
	touch src/RcppExports.cpp

################################################################################
# docker

IMAGE_NAME = cayek/tess3r:latest
CONTAINER_NAME = tess3r

include	docker.mk
