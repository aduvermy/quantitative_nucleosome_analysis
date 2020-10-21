all: TP_experimental_biologists.pdf TP_computational_biologists.pdf ../public/TP_experimental_biologists.html ../public/TP_computational_biologists.html

../public/TP_experimental_biologists.html: TP_experimental_biologists.md
	pandoc -s TP_experimental_biologists.md -o ../public/TP_experimental_biologists.html

../public/TP_computational_biologists.html: TP_computational_biologists.md
	pandoc -s TP_computational_biologists.md -o ../public/TP_computational_biologists.html

TP_experimental_biologists.pdf: TP_experimental_biologists.md
	R -e 'require(rmarkdown); rmarkdown::render("TP_experimental_biologists.md")'

TP_computational_biologists.pdf: TP_computational_biologists.md
	R -e 'require(rmarkdown); rmarkdown::render("TP_computational_biologists.md")'
