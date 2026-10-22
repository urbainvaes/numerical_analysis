# Typst version of the lecture notes: `make` builds ../build/typst/main.pdf
# and one PDF per chapter next to it.

.PHONY: all clean

# Fonts: those built into Typst, plus Font Awesome for icons (Arch package
# otf-font-awesome). figures/ and references.bib are symlinked here.
unexport TYPST_FONT_PATHS
OUT := ../build/typst
SOURCES := $(wildcard *.typ */*.typ) $(wildcard csl/*.csl) ../references.bib

all: $(OUT)/main.pdf

$(OUT)/main.pdf: $(SOURCES)
	mkdir -p $(OUT)
	typst compile main.typ $@
	python3 split.py $@ $(OUT) main.typ

csl/fullcite.csl csl/number.csl: csl/trad-abbrv.csl csl/make_variants.py
	python3 csl/make_variants.py

clean:
	rm -rf $(OUT)
