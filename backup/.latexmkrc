$pdf_mode = 1;
$pdflatex = 'lualatex -synctex=1 --shell-escape --interaction=nonstopmode';
$pdflualatex = 'lualatex -synctex=1 --shell-escape --interaction=nonstopmode';
$out_dir = 'build';
$clean_ext = "synctex.gz synctex.gz(busy) bbl dvi nav run run.xml snm";
$new_viewer_always [0];
$max_repeat = 6;
@default_files = ('main.tex');
