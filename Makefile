.PHONY: all clean

TEX_FILES = main.tex $(wildcard frontmatter/*enpc.tex) $(wildcard mainmatter/*.tex)
PDF_OUTPUTS := $(addprefix build/, $(notdir $(TEX_FILES:.tex=.pdf)))
COMMON_DEPS := main.tex preamble.tex macros.tex $(wildcard mainmatter/**/*.tex)

all : $(PDF_OUTPUTS)

main : build/main.pdf

build/main.pdf: main.tex $(TEX_FILES) $(COMMON_DEPS)
	latexmk -f $<

build/%.pdf: frontmatter/%.tex build/main.pdf $(COMMON_DEPS)
	cp build/main.aux build/$*.aux
	cp build/main.bbl build/$*.bbl
	lualatex \
		--output-directory='build' \
		-shell-escape \
		-jobname=$* \
		"\includeonly{frontmatter/$*}\input{main.tex}"

build/%.pdf: mainmatter/%.tex build/main.pdf $(COMMON_DEPS)
	cp build/main.aux build/$*.aux
	cp build/main.bbl build/$*.bbl
	lualatex \
		--output-directory='build' \
		-shell-escape \
		-jobname=$* \
		"\includeonly{mainmatter/$*}\input{main.tex}"

clean:
	rm -rf build
