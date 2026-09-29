#!/usr/bin/env python3
"""
inspect_gb.py  -  reveal the structure of raghwani_hcv.gb so we can find the
patient/sample grouping. Run it, then paste the output back.

  python inspect_gb.py --gb raghwani_hcv.gb
"""
import argparse, re
from collections import Counter

def split_records(text):
    rec, recs = [], []
    for line in text.splitlines():
        if line.strip() == "//":
            if rec: recs.append("\n".join(rec)); rec = []
        else:
            rec.append(line)
    if rec: recs.append("\n".join(rec))
    return recs

def qual(record, name):
    m = re.search(r'/%s="((?:[^"]|\n)*?)"' % name, record)
    return re.sub(r'\s+', ' ', m.group(1)).strip() if m else ""

ap = argparse.ArgumentParser()
ap.add_argument("--gb", required=True)
args = ap.parse_args()

text = open(args.gb, encoding="utf-8", errors="replace").read()
recs = split_records(text)
print(f"total records: {len(recs)}\n")

print("=" * 70)
print("FIRST FULL RECORD (verbatim):")
print("=" * 70)
print(recs[0][:3000])
print("=" * 70)
print("SECOND FULL RECORD (verbatim):")
print("=" * 70)
print(recs[1][:3000] if len(recs) > 1 else "(only one record)")

# which qualifiers exist
keys = Counter()
for r in recs:
    for k in re.findall(r'/(\w+)=', r):
        keys[k] += 1
print("\n" + "=" * 70)
print("QUALIFIER KEYS (how many records contain each):")
for k, c in keys.most_common():
    print(f"   /{k:<18} {c}")

# DEFINITION lines (top 15 distinct)
defs = Counter(re.sub(r'\s+', ' ', m.group(1)).strip()
               for r in recs for m in [re.search(r'^DEFINITION\s+(.*)', r, re.M)] if m)
print("\nDISTINCT DEFINITION lines:", len(defs), "(top 10)")
for d, c in defs.most_common(10):
    print(f"   [{c}] {d[:90]}")

# values of candidate grouping fields
for field in ("isolate", "strain", "host", "clone", "isolation_source", "note",
              "collection_date", "country", "geo_loc_name"):
    vals = Counter(qual(r, field) for r in recs if qual(r, field))
    if vals:
        print(f"\n/{field}: {len(vals)} distinct values across {sum(vals.values())} records (top 12):")
        for v, c in vals.most_common(12):
            print(f"   [{c}] {v[:60]!r}")

# decompose the isolate by '_' to see which positions vary
isos = [qual(r, "isolate") for r in recs if qual(r, "isolate")]
if isos:
    maxparts = max(len(i.split("_")) for i in isos)
    print(f"\nISOLATE decomposed by '_' (positions, distinct values):")
    for p in range(maxparts):
        vals = Counter(i.split("_")[p] for i in isos if len(i.split("_")) > p)
        sample = ", ".join(list(vals)[:8])
        print(f"   pos {p}: {len(vals)} distinct  e.g. {sample}")
