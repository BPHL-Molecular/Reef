# Reef

Two single-file, browser-based tools for visualising pathogen genomes in public health genomics:

- **Loggerhead** — phylogenetic tree viewer (reads Auspice/Nextstrain JSON or Newick)
- **Calusa** — genetic-distance network viewer (FASTA → `calusa.py` → `network.json`)

**Live:** https://bphl-molecular.github.io/Reef/
· Loggerhead: `…/Reef/loggerhead/` · Calusa: `…/Reef/calusa/calusa.html`

Everything runs in your browser. The sequence file you open is read and rendered on your own
machine — **nothing is uploaded** — and each page works **offline** once saved.

## Repository layout
```
index.html             # hub landing page (this entry point)
loggerhead/index.html  # Loggerhead tool (single file)
calusa/calusa.html     # Calusa viewer (single file)
calusa/calusa.py       # builds network.json from aligned FASTA (Python stdlib)
scripts/               # optional: analysis & pipeline scripts
LICENSE
```

## Use offline
Download or clone the repository and open `index.html` (or `loggerhead/index.html`,
`calusa/calusa.html`) directly in a browser — no server or installation needed.

## Calusa workflow (network)
```
mafft --auto sequences.fasta > aln.fasta
python calusa/calusa.py --input aln.fasta --threshold 0.037 --output ./out
# then open calusa/calusa.html and load out/network.json
```

## Cite
[Authors]. Reef: single-file, browser-based visualisation of pathogen phylogenies and
genetic-distance networks for public health genomics. *Microbial Genomics* (year).
Software: Zenodo DOI [to be minted].

## License
[Choose one, e.g. MIT / BSD-3-Clause / Apache-2.0] — add a LICENSE file.
