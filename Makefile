PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: Rcpp compile

# compiles C++ code and generates the shared library
compile:
	Rscript -e 'pkgbuild::compile_dll(compile_attributes = FALSE, debug=FALSE)'

# generates Rcpp exports
Rcpp:
	Rscript -e "Rcpp::compileAttributes()"
	
roxygen:
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

# Full suite. Tests of the AD workflow and DLL loading only run against an
# installed package (load_all can't reliably exercise the native-pointer/DLL
# lifecycle), so the default target installs with tests and runs test_package.
test:
	R CMD INSTALL --install-tests .
	Rscript -e 'testthat::test_package("$(PACKAGE)")'

# Fast dev loop: compiles and runs via load_all. Skips the AD/DLL-lifecycle
# tests -- use `make test` to exercise those.
test-local: all
	Rscript -e 'testthat::test_local()'

install:
	R CMD INSTALL .

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

clean:
	rm -f src/*.o src/*.so src/*.o.tmp

vignettes:
	Rscript -e "devtools::build_vignettes()"

.PHONY: all clean test test-local roxygen install build check vignettes
