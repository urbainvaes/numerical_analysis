JULIA_SCRIPTS := $(wildcard *.jl)
HTML_OUTPUTS := $(JULIA_SCRIPTS:.jl=.html)

all : $(HTML_OUTPUTS)

%.html : %.jl
	pygmentize -O full -l julia -o $@ $<
