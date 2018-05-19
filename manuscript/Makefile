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

R: $(knitrfile)
	Rscript -e "Stangle('$(knitrfile)')"

clean:
	rm *.aux *.bbl *.log *.blg *.out *.pdf $(texfile)

deepClean: 
	rm -rf cache *.aux *.bbl *.log *.blg *.out *.pdf $(texfile)
