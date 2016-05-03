filestem := clusrank
knitrfile := $(filestem).Rnw
texfile := $(filestem).tex
pdffile := $(filestem).pdf

all: tex pdf

tex : $(knitrfile)
	$(info ************  CONVERTING KNITR TO TEX  ************)
	Rscript -e "require(knitr); knit('$(knitrfile)')"

pdf : tex
	$(info ************  COMPILING PDF  ************)
	pdflatex $(filestem)
	bibtex $(filestem)
	pdflatex $(filestem)
	pdflatex $(filestem)

clean:
	rm *.aux *.bbl *.log *.blg
# watch:
#     $(info ************  WATCHING FOR CHANGES  ************)
#     watchman watch "$(shell pwd)"
#     watchman -- trigger "$(shell pwd)" remake *.Rnw -- make pdf

# unwatch:
#     $(info ************  CANCELLING WATCH  ************)
#     watchman watch-del "$(shell pwd)"
