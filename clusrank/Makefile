objects := $(wildcard R/*.R) DESCRIPTION
dir := $(shell pwd)
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log
yr := $(shell date +"%Y")
dt := $(shell date +"%Y-%m-%d")
cprt := COPYRIGHT


.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

$(tar): $(objects)
	@make -s updateMeta
	Rscript -e "roxygen2::roxygenize(); library(methods); devtools::document();";
	R CMD build $(dir) --compact-vignettes="gs+qpdf"

$(checkLog): $(tar)
	R CMD check --as-cran $(tar)

## update copyright year in R script and date in DESCRIPTION
.PHONY: updateMeta
updateMeta: $(cprt) $(objects)
	@echo "Updating date, version, and copyright year"
	@sed -i "s/Copyright (C) 20[0-9]\{2\}-*[0-9]*/Copyright (C) 2015-$(yr)/" $(cprt)
	@for Rfile in R/*.R; do \
	if ! grep -q 'Copyright (C)' $$Rfile;\
	then cat $(cprt) $$Rfile > tmp;\
	mv tmp $$Rfile;\
	fi;\
	sed -i "s/Copyright (C) 20[0-9]\{2\}-*[0-9]*/Copyright (C) 2015-$(yr)/" $$Rfile;\
	done;
	@sed -i "s/Date: [0-9]\{4\}-[0-9]\{1,2\}-[0-9]\{1,2\}/Date: $(dt)/" DESCRIPTION

.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy src/*.o *.tar.gz *.Rcheck/ .\#*
