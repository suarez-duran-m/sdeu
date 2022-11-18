base=auger_bxl_26oct2022
#base=auger_Wuppertal_21jul2022
#auger_bxl_13jul2022
#bxl_26jan2022
#sdeuMeetOct122021
#auger_bxl
#miniReport
#sdeuMonitoring
#version=_v5
name=${base}
#name=${base}${version}
all:
	pdflatex ${name}.tex

#cuando esta terminado, se copia library.lib, se edita y se comentan los bibtex.

full:
	make clean
	pdflatex ${name}.tex
	bibtex ${name}.aux
	pdflatex ${name}.tex
	pdflatex ${name}.tex
	pdflatex ${name}.tex

clean:
	rm -f ${name}.aux ${name}.log ${name}.bbl ${name}.blg ${name}.out ${name}.toc sections/*.aux ${name}.dvi ${name}.nav ${name}.snm ${name}.vrb
fullclean:
	rm -f *.aux *.log *.gz *.bak *.bbl *.blg *.out *.toc ${name}.pdf
