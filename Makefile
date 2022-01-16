TEX_FILES := main.tex $(notdir frontmatter/syllabus.tex $(wildcard mainmatter/*.tex))
PDF_OUTPUTS := $(addprefix build/, $(TEX_FILES:.tex=.pdf))

all : $(PDF_OUTPUTS)

main : build/main.pdf

build/main.pdf: main.tex
	latexmk $^

build/syllabus.pdf: frontmatter/syllabus.tex build/main.pdf
	lualatex \
		--result=$@ \
		--output-directory='build' \
		-shell-escape \
		-jobname=syllabus \
		"\includeonly{$<}\input{main.tex}"

build/%.pdf: mainmatter/%.tex build/main.pdf
	lualatex \
		--result=$@ \
		--output-directory='build' \
		-shell-escape \
		-jobname=$* \
		"\includeonly{$<}\input{main.tex}"
