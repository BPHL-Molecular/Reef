#!/usr/bin/env python3
"""
concordance.py - quantify Calusa clusters vs a metadata grouping (e.g. patient).
Reads a Calusa network.json and the reef_metadata.tsv; prints a cluster x group
cross-tab, per-cluster purity, and the Adjusted Rand Index. Produces Table 2.

  python concordance.py --calusa calout/network.json --metadata reef_metadata.tsv \
                        --id-col accession --group-col patient
Pure standard library.
"""
import argparse, csv, json
from collections import defaultdict, Counter
from math import comb

ap = argparse.ArgumentParser()
ap.add_argument("--calusa", required=True)
ap.add_argument("--metadata", required=True)
ap.add_argument("--id-col", default="accession")
ap.add_argument("--group-col", default="patient")
args = ap.parse_args()

net = json.load(open(args.calusa, encoding="utf-8"))
nodes = net.get("nodes", [])
node_cluster = {str(n.get("id")): n.get("cluster") for n in nodes}

group = {}
with open(args.metadata, encoding="utf-8-sig", newline="") as fh:
    delim = "\t" if args.metadata.lower().endswith((".tsv", ".txt", ".tab")) else ","
    for row in csv.DictReader(fh, delimiter=delim):
        group[row[args.id_col].strip()] = row[args.group_col].strip()

# pair up only nodes present in both
ids = [i for i in node_cluster if node_cluster[i] is not None and i in group]
ct = defaultdict(Counter)                      # cluster -> {group: n}
for i in ids:
    ct[node_cluster[i]][group[i]] += 1

clusters = sorted(ct, key=lambda c: (str(type(c)), c))
groups = sorted({group[i] for i in ids})

print(f"nodes analysed: {len(ids)}   clusters: {len(clusters)}   {args.group_col}s: {len(groups)}\n")
print("cluster x %s cross-tab:" % args.group_col)
hdr = "cluster".ljust(9) + "".join(g[:8].rjust(9) for g in groups) + "   size  purity dominant"
print(hdr); print("-" * len(hdr))
pure = 0
for c in clusters:
    row = ct[c]; size = sum(row.values())
    dom, dn = row.most_common(1)[0]
    purity = dn / size
    if purity == 1.0: pure += 1
    cells = "".join(str(row.get(g, 0)).rjust(9) for g in groups)
    print(f"{str(c).ljust(9)}{cells}   {str(size).rjust(4)}  {purity:5.0%}  {dom}")
print("-" * len(hdr))
print(f"\npure clusters (subset of one {args.group_col}): {pure}/{len(clusters)}")

# Adjusted Rand Index between cluster labels and group labels
def ari(labels_a, labels_b):
    ct2 = defaultdict(Counter)
    for a, b in zip(labels_a, labels_b): ct2[a][b] += 1
    nij = sum(comb(v, 2) for row in ct2.values() for v in row.values())
    a_sum = sum(comb(sum(row.values()), 2) for row in ct2.values())
    b_tot = Counter()
    for row in ct2.values():
        for b, v in row.items(): b_tot[b] += v
    b_sum = sum(comb(v, 2) for v in b_tot.values())
    n = len(labels_a); tot = comb(n, 2)
    exp = a_sum * b_sum / tot if tot else 0
    maxi = 0.5 * (a_sum + b_sum)
    return (nij - exp) / (maxi - exp) if (maxi - exp) else 1.0

a = ari([node_cluster[i] for i in ids], [group[i] for i in ids])
print(f"Adjusted Rand Index (cluster vs {args.group_col}): {a:.3f}")
print("  ~1.0 = clusters refine groups cleanly; lower = clusters cross group boundaries.")
