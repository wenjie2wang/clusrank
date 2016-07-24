filestem := clusrank
knitrfile := $(filestem).Rnw
texfile := $(filestem).tex
bibfile := $(filestem).bib
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
	rm *.aux *.bbl *.log *.blg *.out *.pdf *.tex
# watch:
#     $(info ************  WATCHING FOR CHANGES  ************)
#     watchman watch "$(shell pwd)"
#     watchman -- trigger "$(shell pwd)" remake *.Rnw -- make pdf

# unwatch:
#     $(info ************  CANCELLING WATCH  ************)
#     watchman watch-del "$(shell pwd)"
