# Reef

**Reef** is a suite of single-file, browser-based tools for visualizing pathogen genomes in public health genomics. It bundles two viewers:

- **Loggerhead** — an interactive phylogenetic **tree** viewer (reads Auspice/Nextstrain JSON or Newick)
- **Calusa** — an interactive genetic-distance **network** viewer (FASTA → `calusa.py` → `network.json`)

Everything runs in your browser. The file you open is read and rendered on your own machine — **nothing is uploaded** — and each page works **offline** once saved. That matters when genomic data is tied to identifiable, reportable conditions and cannot be sent to an outside service.

> **Reef supersedes the previous standalone Calusa project.** Calusa now lives inside Reef, alongside Loggerhead, and all further development happens here. If you were using the old standalone Calusa, switch to this repository — it is the same viewer, now maintained as part of the Reef suite.

**Live:** https://bphl-molecular.github.io/Reef/
&nbsp;·&nbsp; Loggerhead: https://bphl-molecular.github.io/Reef/loggerhead/loggerhead.html
&nbsp;·&nbsp; Calusa: https://bphl-molecular.github.io/Reef/calusa/calusa.html

---

## Loggerhead — phylogenetic tree viewer

Loggerhead renders and explores phylogenies directly in the browser.

- **Reads:** Auspice / Nextstrain v2 JSON (tree + metadata + colorings load together), or plain **Newick** (`.nwk`, `.tree`, `.tre`) from IQ-TREE, RAxML, FastTree, etc.
- **Layouts:** rectangular and radial.
- **Explore:** color tips and branches by metadata, search and highlight tips, and view divergence vs. sampling date for time-resolved trees.
- **Export:** save the current view as PNG or SVG for figures and reports.
- **Demo:** a built-in example tree loads instantly — no input needed.

*Named for the loggerhead sea turtle, which nests along Florida's coasts.*

## Calusa — genetic-distance network viewer

Calusa visualizes pairwise genetic relationships as a network, complementing a tree by showing which samples cluster tightly at a chosen distance cutoff.

- **Distance:** pairwise **GHOST-style corrected Hamming distance** between aligned sequences (corrected for comparable positions rather than a raw mismatch count).
- **Reads:** a `network.json` produced by `calusa/calusa.py` from an aligned FASTA.
- **Shows:** a force-directed graph with distance-threshold clusters and link weighting.
- **Threshold:** adjust the cutoff to watch clusters merge and split — useful for choosing a defensible threshold.
- **Export:** save the network as PNG or SVG.
- **Demo:** a built-in example network loads instantly — no input needed.

> Calusa is a **genetic-distance** viewer: distance-based clusters describe sequence similarity, which is evidence to interpret alongside epidemiological data — not a direct claim of who infected whom.

*Named for Calusa Beach on Bahia Honda Key, Florida.*

---

## Repository layout

```
index.html               # hub landing page (entry point)
loggerhead/
  loggerhead.html        # Loggerhead — tree viewer (single file)
calusa/
  calusa.html            # Calusa — network viewer (single file)
  calusa.py              # builds network.json from aligned FASTA (Python stdlib only)
LICENSE                  # MIT
README.md
```

## Use offline

Download or clone the repository and open `index.html` (or `loggerhead/loggerhead.html`, `calusa/calusa.html`) directly in a browser — no server or installation needed.

```bash
git clone https://github.com/BPHL-Molecular/Reef.git
cd Reef
# then open index.html in your browser
```

## Calusa workflow (building a network)

`calusa.py` needs only the Python standard library.

```bash
# 1. align your sequences
mafft --auto sequences.fasta > aln.fasta

# 2. build the network (default threshold 0.037)
python calusa/calusa.py --input aln.fasta --threshold 0.037 --output ./out

# 3. open calusa/calusa.html and load out/network.json
```

`calusa.py` options: `--input/-i` (aligned FASTA), `--threshold/-t` (distance cutoff, default `0.037`), `--output/-o` (output directory, default `.`), and `--create-sample` to write an example FASTA.



## License

Released under the **MIT License** — see [LICENSE](LICENSE). You may use, modify, and redistribute Reef, including commercially, provided the copyright and license notice are retained.
