
clean:
	rm -f *.log *.aux *.bbl *.blg *.pdf *.gz *.backup *.bak *.toc *.idx *.ilg *.ind

check: 
	grep -e 'TODO|FIXME' *.tex
	grep -Eo '(\b.+) \1\b' *.tex 

spelling: 
	for Datei in *.tex; do ispell -t $$Datei; done

quick: 
	pdflatex -synctex=1 -interaction=nonstopmode main.matrices.tex

main: 
	pdflatex -synctex=1 -interaction=nonstopmode main.matrices.tex 
	bibtex main.matrices.aux
	pdflatex -synctex=1 -interaction=nonstopmode main.matrices.tex 
	pdflatex -synctex=1 -interaction=nonstopmode main.matrices.tex

