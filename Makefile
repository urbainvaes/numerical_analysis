all : build/main.pdf build/roundoff.pdf

build/main.pdf: main.tex
	lualatex \
		--result=$@ \
		--output-directory='build' \
		--shell-escape \
		main.tex


build/roundoff.pdf: mainmatter/roundoff.tex
	lualatex \
		--result=$@ \
		--output-directory='build' \
		-shell-escape \
		-jobname=roundoff \
		"\includeonly{mainmatter/roundoff.tex}\input{main.tex}"
