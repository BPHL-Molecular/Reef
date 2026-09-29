#!/usr/bin/env python3
r"""
prep_reef_inputs.py  (v2)
=========================
Turn the downloaded GenBank file into Reef inputs. Works for the Raghwani 2016
HCV data, whose /isolate is  {variant}_{patient}_{timepoint}_{offset}
(e.g. "1_russo_3_27"): patient = field index 1.

Patient is taken from --patient-regex (capture group 1) applied to /isolate.
  For the Raghwani data:  --patient-regex '^\d+_([^_]+)'   ->  patient = russo, ricci, ...

TWO MODES
  (default) grouped : one node/tip per PATIENT (8 nodes). calusa headers {patient}_seq_k;
                      tree = one representative sequence per patient.
  --ungrouped       : one node/tip per VARIANT (recommended for this intra-host dataset).
                      calusa headers = accession; tree = every variant (header = accession).
                      Expect ~8 clusters (one per patient) at threshold 0.037.
  --max-per-patient N : cap variants per patient (keeps figures legible; e.g. 30).

  python prep_reef_inputs.py --gb raghwani_hcv.gb --inspect
  python prep_reef_inputs.py --gb raghwani_hcv.gb --patient-regex '^\d+_([^_]+)' --ungrouped --max-per-patient 30

Pure standard library.
"""
import argparse, re, sys
from collections import defaultdict, Counter
import datetime as dt

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

def get_accession(record):
    m = re.search(r'^ACCESSION\s+(\S+)', record, re.M)
    return m.group(1) if m else ""

def get_sequence(record):
    m = re.search(r'\nORIGIN(.*)$', record, re.S)
    return re.sub(r'[^A-Za-z]', '', m.group(1)).upper() if m else ""

def norm_date(s):
    s = s.strip()
    if not s: return ""
    for fmt, out in (("%d-%b-%Y", "%Y-%m-%d"), ("%b-%Y", "%Y-%m"), ("%Y-%m-%d", "%Y-%m-%d"),
                     ("%Y-%m", "%Y-%m"), ("%Y", "%Y")):
        try: return dt.datetime.strptime(s, fmt).strftime(out)
        except ValueError: continue
    m = re.search(r'(\d{4})', s)
    return m.group(1) if m else ""

def get_genotype(record):
    g = qual(record, "genotype")
    if g: return g
    m = re.search(r'(?:genotype|subtype)[:\s]+([0-9][a-zA-Z]?)', qual(record, "note"), re.I)
    if m: return m.group(1)
    m2 = re.search(r'subtype\s+([0-9][a-zA-Z]?)', qual(record, "organism"), re.I)
    return m2.group(1) if m2 else ""

def patient_of(isolate, accession, regex):
    if not isolate: return accession
    if regex:
        m = regex.search(isolate)
        if m: return m.group(1) if m.groups() else m.group(0)
        return isolate
    return isolate

def consensus_pick(seqs):
    counts = Counter(seqs)
    return sorted(counts.items(), key=lambda kv: (-kv[1], -len(kv[0])))[0][0]

def mode(values):
    vals = [v for v in values if v]
    return Counter(vals).most_common(1)[0][0] if vals else ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gb", required=True)
    ap.add_argument("--patient-regex", default=None, help="group 1 = patient ID, applied to /isolate")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--ungrouped", action="store_true", help="one node/tip per variant (node id = accession)")
    ap.add_argument("--max-per-patient", type=int, default=0, help="cap variants per patient (0 = no cap)")
    ap.add_argument("--out-calusa", default="calusa_input.fasta")
    ap.add_argument("--out-tree", default="tree_input.fasta")
    ap.add_argument("--out-meta", default="reef_metadata.tsv")
    args = ap.parse_args()

    regex = re.compile(args.patient_regex) if args.patient_regex else None
    records = split_records(open(args.gb, encoding="utf-8", errors="replace").read())
    if not records: sys.exit("ERROR: no records parsed.")

    rows = []
    for rec in records:
        seq = get_sequence(rec)
        if not seq: continue
        iso = qual(rec, "isolate") or qual(rec, "strain")
        parts = iso.split("_")
        rows.append({
            "acc": get_accession(rec), "isolate": iso, "seq": seq,
            "patient": patient_of(iso, get_accession(rec), regex),
            "timepoint": parts[2] if len(parts) > 2 else "",
            "genotype": get_genotype(rec),
            "country": qual(rec, "country") or qual(rec, "geo_loc_name"),
            "date": norm_date(qual(rec, "collection_date")),
        })

    by_patient = defaultdict(list)
    for r in rows: by_patient[r["patient"]].append(r)

    examples = list({r["isolate"]: r["patient"] for r in rows}.items())[:12]
    sizes = sorted((len(v) for v in by_patient.values()), reverse=True)
    print(f"records with sequence : {len(rows)}")
    print(f"distinct patients     : {len(by_patient)}  -> {sorted(by_patient)}")
    print("example isolate -> patient:")
    for iso, pat in examples: print(f"   {iso!r:24s} -> {pat!r}")
    print(f"variants per patient  : max {sizes[0]}, median {sizes[len(sizes)//2]}, total {sum(sizes)}")
    if args.inspect:
        print("\n(--inspect) no files written.")
        return
    if len(by_patient) == len(rows):
        print("\nWARNING: no grouping (every sequence is its own patient). Set --patient-regex.")

    # optional cap per patient
    if args.max_per_patient > 0:
        for p in by_patient: by_patient[p] = by_patient[p][:args.max_per_patient]
    kept = [r for recs in by_patient.values() for r in recs]

    if args.ungrouped:
        with open(args.out_calusa, "w") as f:
            for r in kept: f.write(f">{r['acc']}\n{r['seq']}\n")
        with open(args.out_tree, "w") as f:
            for r in kept: f.write(f">{r['acc']}\n{r['seq']}\n")
        with open(args.out_meta, "w") as f:
            f.write("accession\tpatient\ttimepoint\tgenotype\tcountry\tcollection_date\n")
            for r in kept:
                f.write("\t".join([r["acc"], r["patient"], r["timepoint"],
                                   r["genotype"], r["country"], r["date"]]) + "\n")
        idcol, cols = "accession", "patient,genotype,country"
    else:
        with open(args.out_calusa, "w") as f:
            for pat, recs in by_patient.items():
                for k, r in enumerate(recs, 1): f.write(f">{pat}_seq_{k}\n{r['seq']}\n")
        with open(args.out_tree, "w") as f:
            for pat, recs in by_patient.items():
                f.write(f">{pat}\n{consensus_pick([r['seq'] for r in recs])}\n")
        with open(args.out_meta, "w") as f:
            f.write("patient\tgenotype\tcountry\tcollection_date\tn_variants\n")
            for pat, recs in by_patient.items():
                f.write("\t".join([pat, mode([r['genotype'] for r in recs]),
                                   mode([r['country'] for r in recs]),
                                   mode([r['date'] for r in recs]), str(len(recs))]) + "\n")
        idcol, cols = "patient", "genotype,country"

    print(f"\nwrote {args.out_calusa}, {args.out_tree}, {args.out_meta}  (kept {len(kept)} sequences)")
    print("next (after mafft-aligning both FASTAs and running calusa.py):")
    print(f"   reef_newick_to_auspice.py --tree tree.nwk --metadata {args.out_meta} "
          f"--id-col {idcol} --date-col collection_date --cols {cols} "
          f"--calusa calout/network.json --output reef_auspice.json")

if __name__ == "__main__":
    main()
