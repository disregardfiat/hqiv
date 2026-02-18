# HQIV Paper — LaTeX

- **main.tex** — Full manuscript (abstract, §§1–8, acknowledgments, bibliography).
- **refs.bib** — McCulloch QI and related references.

Build:

```bash
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

Or upload `main.tex` and `refs.bib` to Overleaf; use **natbib** and a style such as `mnras` (or `apj`, `aasjournal` as desired). arXiv: upload the same two files; biblatex users can switch to `biblatex` + `biber` if preferred.
