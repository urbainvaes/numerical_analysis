"""Extract one PDF per part of the book, using the <part> markers of main.typ.

Usage: python3 split.py BOOK.pdf OUTPUT_DIR TYPST_QUERY_ARGS...
"""

import json
import subprocess
import sys

book, out_dir, query_args = sys.argv[1], sys.argv[2], sys.argv[3:]

query = subprocess.run(
    ["typst", "query", *query_args, "<part>", "--field", "value"],
    check=True, capture_output=True, text=True,
)
parts = json.loads(query.stdout)
n_pages = int(subprocess.run(
    ["qpdf", "--show-npages", book], check=True, capture_output=True, text=True,
).stdout)

for i, part in enumerate(parts):
    first = part["page"]
    last = parts[i + 1]["page"] - 1 if i + 1 < len(parts) else n_pages
    if part["output"]:
        output = f"{out_dir}/{part['name']}.pdf"
        subprocess.run(["qpdf", book, "--pages", ".", f"{first}-{last}", "--", output], check=True)
        print(f"{output}: pages {first}-{last}")
