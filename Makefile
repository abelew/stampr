VERSION=2019.12
export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean roxygen reference check build test

build:
	@echo "Performing build with R CMD build errRt"
	R CMD build .

check: roxygen
	@echo "Performing check with R CMD check errRt"
	rm -rf ./..Rcheck 2>/dev/null 1>&2
	export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check . --no-build-vignettes
	@rm -rf ./..Rcheck 2>/dev/null 1>&2

clean:
	@echo "Cleaning up"
	rm -rf errRt
	rm -rf ./..Rcheck &
	rm -rf errRt.Rcheck/
	rm -f errRt_${VERSION}.tar.gz

clean_vignette:
	rm -f vignettes/*.rda vignettes/*.map vignettes/*.Rdata

deps:
	@echo "Invoking devtools::install_dev_deps()"
	R -e "all = as.data.frame(devtools::dev_package_deps('.', dependencies=TRUE)); needed = all[['diff']] < 0; needed = all[needed, 'package']; for (t in needed) { BiocManager::install(t) }"

document: roxygen vignette reference

install:
	@echo "Performing R CMD INSTALL errRt."
	R CMD INSTALL --install-tests .


reference:
	@echo "Generating reference manual with R CMD Rd2pdf"
	@mkdir -p inst/doc
	R CMD Rd2pdf . -o inst/reference/reference.pdf --no-preview

roxygen:
	@echo "Generating documentation with devtools::document()"
	R -e "suppressPackageStartupMessages(devtools::document())"

test: install
	@echo "Running run_tests.R"
	tests/testthat.R

vignette:
	@mkdir -p doc
	@echo "Building vignettes with devtools::build_vignettes()"
	R -e "devtools::build_vignettes()"
	mv doc inst/doc
	cp inst/reference/* inst/doc

vt:	clean_vignette vignette reference install
