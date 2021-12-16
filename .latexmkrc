$pdf_mode = 1;
$pdflatex = 'lualatex -synctex=1 --shell-escape';
$pdflualatex = 'lualatex -synctex=1 --shell-escape';
$out_dir = 'build';
$clean_ext = "synctex.gz synctex.gz(busy) bbl dvi nav run run.xml snm";
$new_viewer_always [0];
@default_files = ('main.tex');
