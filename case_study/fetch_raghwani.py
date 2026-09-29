#!/usr/bin/env python3
"""
fetch_raghwani.py
=================
Download the Raghwani et al. 2016 (PLOS Pathogens) HCV sequences from NCBI.
Accession range: KX111382-KX113360 (1,979 records).

Saves both:
  - raghwani_hcv.fasta   (sequences, for the tree + Calusa)
  - raghwani_hcv.gb      (GenBank flatfile, for metadata: isolate/date/country/genotype)

RUN THIS ON A MACHINE WITH INTERNET (NCBI is not reachable from the Reef sandbox).
NCBI etiquette: supply your email; an API key raises the rate limit to 10 req/s.

  python fetch_raghwani.py --email you@dept.org
  python fetch_raghwani.py --email you@dept.org --api-key YOUR_NCBI_KEY

Pure standard library.
"""
import argparse, sys, time, urllib.parse, urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
START, END, PREFIX = 111382, 113360, "KX"   # KX111382 .. KX113360

def accessions():
    return [f"{PREFIX}{n}" for n in range(START, END + 1)]

def fetch_batch(ids, rettype, email, api_key, retries=4):
    params = {"db": "nuccore", "id": ",".join(ids),
              "rettype": rettype, "retmode": "text", "tool": "reef", "email": email}
    if api_key:
        params["api_key"] = api_key
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(EFETCH, data=data, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            wait = 2 * (attempt + 1)
            sys.stderr.write(f"  retry {attempt+1}/{retries} after error: {e} (wait {wait}s)\n")
            time.sleep(wait)
    raise RuntimeError("batch failed after retries")

def download(rettype, outpath, ids, email, api_key, batch, pause):
    n = len(ids)
    with open(outpath, "w", encoding="utf-8") as out:
        for i in range(0, n, batch):
            chunk = ids[i:i + batch]
            txt = fetch_batch(chunk, rettype, email, api_key)
            out.write(txt)
            if not txt.endswith("\n"):
                out.write("\n")
            print(f"  {rettype}: {min(i+batch, n)}/{n}")
            time.sleep(pause)
    print(f"wrote {outpath}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="your email (NCBI requirement)")
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--fasta-out", default="raghwani_hcv.fasta")
    ap.add_argument("--gb-out", default="raghwani_hcv.gb")
    ap.add_argument("--batch", type=int, default=200)
    args = ap.parse_args()

    pause = 0.12 if args.api_key else 0.4   # ~10/s with key, ~3/s without
    ids = accessions()
    print(f"Downloading {len(ids)} records ({ids[0]}..{ids[-1]})")
    download("fasta", args.fasta_out, ids, args.email, args.api_key, args.batch, pause)
    download("gb",    args.gb_out,    ids, args.email, args.api_key, args.batch, pause)
    print("done.")

if __name__ == "__main__":
    main()
