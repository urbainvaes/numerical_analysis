JULIA_SCRIPTS := $(wildcard notebook*.jl)
IPYNB_OUTPUTS := $(JULIA_SCRIPTS:.jl=.ipynb)
HTML_OUTPUTS := $(JULIA_SCRIPTS:.jl=.html)
CSS := .jupyter/custom/custom.css

all : $(HTML_OUTPUTS) $(IPYNB_OUTPUTS)

%.ipynb : %.jl $(CSS)
	~/.local/bin/jupytext --to 'notebook' --set-kernel 'julia-1.5' $<
	JUPYTER_CONFIG_DIR=.jupyter jupyter nbconvert --execute --to 'notebook' --inplace $@

%.html : %.ipynb $(CSS)
	JUPYTER_CONFIG_DIR=.jupyter jupyter nbconvert --to html $<

clean :
	rm $(IPYNB_OUTPUTS) $(HTML_OUTPUTS)
