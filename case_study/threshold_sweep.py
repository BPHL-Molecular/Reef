#!/usr/bin/env python3
"""
threshold_sweep.py - Reef Fig S2 (threshold sensitivity)
Runs calusa.py over a range of distance thresholds on the SAME aligned FASTA,
records cluster / link counts, writes a CSV and a 'clusters vs threshold' plot.

  python threshold_sweep.py --input aln.fasta --start 0.02 --stop 0.06 --step 0.005
  # or explicit:  --thresholds 0.0,0.01,0.02,0.037,0.05

Needs calusa.py reachable (default: ./calusa.py) and Python. Plot needs matplotlib
(falls back to CSV-only if absent). Pure standard library otherwise.
"""
import argparse, csv, json, os, subprocess, sys

def cluster_stats(network_json):
    d = json.load(open(network_json, encoding="utf-8"))
    nodes = d.get("nodes", [])
    clusters = {n.get("cluster") for n in nodes if n.get("cluster") is not None}
    linked = sum(1 for n in nodes if n.get("cluster") is not None)
    return len(nodes), linked, len(clusters), len(d.get("links", []))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="aligned FASTA (same one used for the network)")
    ap.add_argument("--calusa-script", default="calusa.py")
    ap.add_argument("--thresholds", help="comma-separated thresholds (overrides start/stop/step)")
    ap.add_argument("--start", type=float, default=0.02)
    ap.add_argument("--stop", type=float, default=0.06)
    ap.add_argument("--step", type=float, default=0.005)
    ap.add_argument("--workdir", default="sweep_out")
    ap.add_argument("--out-csv", default="threshold_sweep.csv")
    ap.add_argument("--out-plot", default="fig_s2_threshold_sweep.png")
    args = ap.parse_args()

    if args.thresholds:
        ts = [round(float(x), 6) for x in args.thresholds.split(",")]
    else:
        ts, t = [], args.start
        while t <= args.stop + 1e-9:
            ts.append(round(t, 6)); t += args.step
    os.makedirs(args.workdir, exist_ok=True)

    rows = []
    for t in ts:
        outdir = os.path.join(args.workdir, f"t_{t}")
        cmd = [sys.executable, args.calusa_script, "--input", args.input,
               "--threshold", str(t), "--output", outdir]
        print(f"threshold {t} ...", flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.stderr.write(r.stdout + r.stderr + "\n")
            sys.exit(f"calusa.py failed at threshold {t}")
        nj = os.path.join(outdir, "network.json")
        n, linked, nclust, nlinks = cluster_stats(nj)
        rows.append({"threshold": t, "n_samples": n, "n_linked": linked,
                     "n_clusters": nclust, "n_links": nlinks})

    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["threshold", "n_samples", "n_linked", "n_clusters", "n_links"])
        w.writeheader(); w.writerows(rows)

    print("\nthreshold  clusters  linked  links")
    for r in rows:
        print(f"  {r['threshold']:<8} {r['n_clusters']:>7}  {r['n_linked']:>6}  {r['n_links']:>6}")
    print(f"\nwrote {args.out_csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xs = [r["threshold"] for r in rows]
        plt.figure(figsize=(6, 4), dpi=300)
        plt.plot(xs, [r["n_clusters"] for r in rows], "o-", color="#1e40af")
        plt.xlabel("Genetic-distance threshold"); plt.ylabel("Number of clusters")
        plt.title("Calusa cluster count vs threshold")
        plt.grid(True, alpha=0.3)
        plt.tight_layout(); plt.savefig(args.out_plot, dpi=300)
        print(f"wrote {args.out_plot}")
    except Exception as e:
        print(f"(plot skipped: {e}) - CSV is available for plotting elsewhere.")

if __name__ == "__main__":
    main()
