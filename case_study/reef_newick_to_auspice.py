#!/usr/bin/env python3
"""
reef_newick_to_auspice.py
=========================

Turn a bare Newick tree into an Auspice v2 JSON that Loggerhead can colour by
metadata (sub-genotype, date, county, ...) and, optionally, by Calusa cluster.

WHY: Newick carries only topology + branch lengths + tip names. Loggerhead builds
its "Color by" menu from an Auspice JSON's `meta.colorings` and each node's
`node_attrs`. A bare Newick therefore offers only "Uniform". This script attaches
the metadata so the colour options appear.

USAGE
-----
  python reef_newick_to_auspice.py \
      --tree hav.nwk \
      --metadata metadata.tsv \
      --id-col strain \
      --date-col date \
      --cols subgenotype,county \
      --calusa network.json \
      --output hav_auspice.json

  - --metadata : CSV/TSV with one row per sequence. Must contain the ID column
                 (matching the Newick tip names) plus your attribute columns.
  - --id-col   : the column holding tip names (default: tries strain/name/id/first).
  - --date-col : a date column -> converted to num_date (decimal year, continuous).
  - --cols     : comma-separated categorical columns to include (default: all others).
  - --calusa   : a Calusa network.json -> adds a `cluster` colouring (for Fig 6).
  - At least one of --metadata / --calusa is required.

Load the resulting *.json in Loggerhead (it accepts .json). Pure standard library.
"""

import argparse, csv, json, sys, datetime as dt

# ---------------------------------------------------------------- Newick parser
def parse_newick(text):
    # strip [comments] that are not inside single quotes
    out, q, i, n = [], False, 0, len(text)
    while i < n:
        c = text[i]
        if c == "'":
            q = not q; out.append(c)
        elif c == '[' and not q:
            depth = 1; i += 1
            while i < n and depth > 0:
                if text[i] == '[': depth += 1
                elif text[i] == ']': depth -= 1
                i += 1
            continue
        else:
            out.append(c)
        i += 1
    s = ''.join(out).strip()
    if s.endswith(';'):
        s = s[:-1]
    pos = 0

    def skip_ws():
        nonlocal pos
        while pos < len(s) and s[pos].isspace():
            pos += 1

    def read_name():
        nonlocal pos
        skip_ws()
        if pos < len(s) and s[pos] == "'":
            pos += 1; name = ''
            while pos < len(s):
                if s[pos] == "'":
                    if pos + 1 < len(s) and s[pos + 1] == "'":
                        name += "'"; pos += 2; continue
                    pos += 1; break
                name += s[pos]; pos += 1
            return name
        name = ''
        while pos < len(s) and s[pos] not in ':,()[;':
            name += s[pos]; pos += 1
        return name.strip()

    def read_len():
        nonlocal pos
        skip_ws()
        if pos < len(s) and s[pos] == ':':
            pos += 1; num = ''
            while pos < len(s) and (s[pos].isdigit() or s[pos] in '+-.eE'):
                num += s[pos]; pos += 1
            try:
                return float(num)
            except ValueError:
                return 0.0
        return 0.0

    def parse_clade():
        nonlocal pos
        node = {"name": "", "bl": 0.0, "children": []}
        skip_ws()
        if pos < len(s) and s[pos] == '(':
            pos += 1
            while True:
                node["children"].append(parse_clade())
                skip_ws()
                if pos < len(s) and s[pos] == ',':
                    pos += 1; continue
                break
            skip_ws()
            if pos < len(s) and s[pos] == ')':
                pos += 1
        node["name"] = read_name()
        node["bl"] = read_len()
        return node

    root = parse_clade()
    root["bl"] = 0.0
    return root

# ---------------------------------------------------------------- helpers
def decimal_year(value):
    s = str(value).strip()
    if not s:
        return None
    if len(s) == 4 and s.isdigit():
        return float(s) + 0.5
    d = None
    for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y", "%d/%m/%Y", "%Y-%m"):
        try:
            d = dt.datetime.strptime(s, fmt).date(); break
        except ValueError:
            continue
    if d is None:
        return None
    y = d.year
    start = dt.date(y, 1, 1).toordinal()
    end = dt.date(y + 1, 1, 1).toordinal()
    return round(y + (d.toordinal() + 0.5 - start) / (end - start), 4)

def load_table(path):
    with open(path, newline='', encoding='utf-8-sig') as fh:
        sample = fh.read(4096); fh.seek(0)
        if path.lower().endswith(('.tsv', '.tab', '.txt')):
            delim = '\t'
        elif path.lower().endswith('.csv'):
            delim = ','
        else:
            delim = '\t' if sample.count('\t') >= sample.count(',') else ','
        return list(csv.DictReader(fh, delimiter=delim))

def pick_id_col(fieldnames, requested):
    if requested:
        if requested not in fieldnames:
            sys.exit(f"ERROR: --id-col '{requested}' not in metadata columns: {fieldnames}")
        return requested
    for cand in ('strain', 'name', 'id', 'sample', 'accession', 'taxon'):
        for f in fieldnames:
            if f.lower() == cand:
                return f
    return fieldnames[0]

# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="Newick + metadata -> Auspice v2 JSON for Loggerhead")
    ap.add_argument('--tree', required=True)
    ap.add_argument('--metadata')
    ap.add_argument('--id-col')
    ap.add_argument('--date-col')
    ap.add_argument('--cols', help="comma-separated categorical columns (default: all non-id/date)")
    ap.add_argument('--calusa', help="Calusa network.json -> adds a cluster colouring")
    ap.add_argument('--title', default="Reef phylogeny")
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if not args.metadata and not args.calusa:
        sys.exit("ERROR: supply --metadata and/or --calusa (otherwise nothing to colour by).")

    root = parse_newick(open(args.tree, encoding='utf-8').read())

    # cumulative divergence
    def assign_div(node, parent=0.0):
        node["_div"] = parent + node.get("bl", 0.0)
        for c in node["children"]:
            assign_div(c, node["_div"])
    assign_div(root, 0.0)
    root["_div"] = 0.0

    # metadata
    meta_map, cat_cols, has_date = {}, [], False
    if args.metadata:
        rows = load_table(args.metadata)
        if not rows:
            sys.exit("ERROR: metadata file is empty.")
        fields = list(rows[0].keys())
        idc = pick_id_col(fields, args.id_col)
        datec = args.date_col
        if datec and datec not in fields:
            sys.exit(f"ERROR: --date-col '{datec}' not in metadata columns: {fields}")
        if args.cols:
            cat_cols = [c.strip() for c in args.cols.split(',') if c.strip()]
            for c in cat_cols:
                if c not in fields:
                    sys.exit(f"ERROR: column '{c}' not in metadata columns: {fields}")
        else:
            cat_cols = [f for f in fields if f != idc and f != datec]
        for r in rows:
            name = (r.get(idc) or '').strip()
            if not name:
                continue
            attrs = {}
            for c in cat_cols:
                v = (r.get(c) or '').strip()
                if v:
                    attrs[c] = v
            if datec:
                dy = decimal_year(r.get(datec, ''))
                if dy is not None:
                    attrs['num_date'] = dy; has_date = True
            meta_map[name] = attrs

    # calusa clusters
    cluster_map = {}
    if args.calusa:
        net = json.load(open(args.calusa, encoding='utf-8'))
        nodes = net.get('nodes', net if isinstance(net, list) else [])
        for nd in nodes:
            nid = nd.get('id') or nd.get('name')
            cl = nd.get('cluster')
            if nid is not None and cl is not None:
                cluster_map[str(nid)] = f"Cluster {cl}"

    # build auspice tree + match stats
    stats = {"tips": 0, "meta_hit": 0, "cluster_hit": 0, "miss": []}

    def to_auspice(node):
        na = {"div": round(node["_div"], 8)}
        name = (node["name"] or "").strip()
        is_tip = not node["children"]
        if is_tip:
            stats["tips"] += 1
            md = meta_map.get(name)
            if md:
                stats["meta_hit"] += 1
                for k, v in md.items():
                    na[k] = {"value": v}
            elif args.metadata and len(stats["miss"]) < 8:
                stats["miss"].append(name)
            if name in cluster_map:
                stats["cluster_hit"] += 1
                na["cluster"] = {"value": cluster_map[name]}
        out = {"name": name if name else None, "node_attrs": na}
        if node["children"]:
            out["children"] = [to_auspice(c) for c in node["children"]]
        return out

    tree = to_auspice(root)

    # colorings
    colorings = []
    for c in cat_cols:
        colorings.append({"key": c, "title": c.replace('_', ' ').title(), "type": "categorical"})
    if has_date:
        colorings.append({"key": "num_date", "title": "Date", "type": "continuous"})
    if cluster_map:
        colorings.append({"key": "cluster", "title": "Calusa cluster", "type": "categorical"})

    default_color = colorings[0]["key"] if colorings else ""
    doc = {
        "version": "v2",
        "meta": {
            "title": args.title,
            "updated": dt.date.today().isoformat(),
            "colorings": colorings,
            "display_defaults": {"color_by": default_color},
        },
        "tree": tree,
    }
    json.dump(doc, open(args.output, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

    # report
    print(f"tips in tree           : {stats['tips']}")
    if args.metadata:
        print(f"tips matched metadata  : {stats['meta_hit']}")
        if stats["miss"]:
            print(f"  unmatched (first few): {', '.join(stats['miss'])}")
            print("  -> check that the metadata ID column matches the Newick tip names exactly.")
    if args.calusa:
        print(f"tips matched a cluster : {stats['cluster_hit']}")
    print(f"colour-by options added: {[c['key'] for c in colorings]}")
    print(f"written                : {args.output}")

if __name__ == '__main__':
    main()
