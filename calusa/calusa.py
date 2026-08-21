#!/usr/bin/env python3
"""
Real Data Processor
Generate threshold_links.csv, cluster_summary.csv, and network.json
from real FASTA sequence files.

Workflow:
1. Read FASTA format sequence files
2. Calculate Hamming distance between all sample pairs
3. Generate transmission network links below threshold
4. Identify transmission clusters
5. Export JSON file for D3.js visualization in transmission-network.html

Usage:
    python real_data_processor.py --input sequences.fasta
    python real_data_processor.py --input sequences.fasta --threshold 0.05 --output ./results
    python real_data_processor.py --create-sample
"""

import json
import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict


def parse_fasta(fasta_file):
    """
    Parse FASTA file.

    Expected format (no frequency info in header):
        >SampleID_seqN
        ATCGATCG...

    Each unique header prefix (everything before the last '_seqN' segment)
    is treated as the sample ID.  All sequences belonging to the same sample
    are grouped together.

    Returns: {sample_id: [sequence, ...]}
    """
    sequences = defaultdict(list)
    current_sample = None
    current_seq = []

    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # Save previous sequence
                    if current_sample and current_seq:
                        sequences[current_sample].append(''.join(current_seq))
                        current_seq = []

                    # Parse header: strip trailing _seqN / _hapN part
                    # e.g. ">Sample_001_seq_1"  -> sample_id = "Sample_001"
                    #      ">Pat123_hap2"        -> sample_id = "Pat123"
                    header = line[1:]
                    parts = header.split('_')

                    if len(parts) > 1 and (parts[-1].isdigit() or
                                           parts[-1].lower().startswith(('seq', 'hap'))):
                        stripped = parts[:]
                        if stripped[-1].isdigit():
                            stripped.pop()
                        if stripped and stripped[-1].lower() in ('seq', 'hap', 'sequence', 'haplotype'):
                            stripped.pop()
                        sample_id = '_'.join(stripped) if stripped else header
                    else:
                        sample_id = header

                    current_sample = sample_id
                else:
                    current_seq.append(line.upper())

            # Save last sequence
            if current_sample and current_seq:
                sequences[current_sample].append(''.join(current_seq))

    except FileNotFoundError:
        print(f"❌ Error: File not found: {fasta_file}")
        return None

    return dict(sequences)


def calculate_hamming_distance(seq1, seq2):
    """
    Corrected Hamming distance (per GHOST manual §5.9.4.c):
    - Only ATCG positions counted
    - Gaps / ambiguous bases skipped
    """
    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]

    differences = 0
    valid_positions = 0

    for c1, c2 in zip(seq1, seq2):
        if c1 in 'ATCG' and c2 in 'ATCG':
            valid_positions += 1
            if c1 != c2:
                differences += 1

    if valid_positions == 0:
        return 1.0

    return differences / valid_positions


def calculate_pairwise_distances(sequences_dict, threshold=0.037):
    """
    Calculate minimum pairwise Hamming distance between all sample pairs.

    Per GHOST manual §5.10.2:
    - Compare every sequence in sample A against every sequence in sample B
    - Take the minimum distance as the distance between the two samples
    - Keep only pairs whose minimum distance is <= threshold
    - For HCV, threshold=0.037
    - For HAV-VPB, threshold=0

    Returns: DataFrame with columns
        source, target, minimum_distance, num_seqs_source, num_seqs_target
    """
    links = []
    sample_ids = list(sequences_dict.keys())

    for i in range(len(sample_ids)):
        for j in range(i + 1, len(sample_ids)):
            s1 = sample_ids[i]
            s2 = sample_ids[j]

            seqs1 = sequences_dict[s1]
            seqs2 = sequences_dict[s2]

            min_dist = 1.0
            for seq1 in seqs1:
                for seq2 in seqs2:
                    d = calculate_hamming_distance(seq1, seq2)
                    if d < min_dist:
                        min_dist = d
                    if min_dist == 0.0:
                        break
                if min_dist == 0.0:
                    break

            if min_dist <= threshold:
                links.append({
                    'source': s1,
                    'target': s2,
                    'minimum_distance': round(min_dist, 6),
                    'num_seqs_source': len(seqs1),
                    'num_seqs_target': len(seqs2),
                })

    return pd.DataFrame(links)


def identify_clusters_from_links(links_df, all_sample_ids):
    """
    Identify transmission clusters using depth-first search on the link graph.

    Parameters:
        links_df       : DataFrame with 'source' and 'target' columns
        all_sample_ids : list of all sample IDs (including unlinked ones)

    Returns:
        dict mapping sample_id -> cluster_id (int), unlinked samples absent
    """
    adjacency = {sid: set() for sid in all_sample_ids}
    for _, row in links_df.iterrows():
        adjacency[row['source']].add(row['target'])
        adjacency[row['target']].add(row['source'])

    visited = set()
    cluster_id = 0
    sample_to_cluster = {}

    def dfs(node, cid):
        visited.add(node)
        sample_to_cluster[node] = cid
        for neighbor in adjacency[node]:
            if neighbor not in visited:
                dfs(neighbor, cid)

    for sid in all_sample_ids:
        if sid not in visited and len(adjacency[sid]) > 0:
            dfs(sid, cluster_id)
            cluster_id += 1

    return sample_to_cluster


def export_network_json(sequences_dict, links_df, sample_to_cluster, threshold):
    """
    Build a JSON structure compatible with transmission-network.html.

    JSON schema
    -----------
    {
      "metadata": {
        "threshold": <float>,
        "total_samples": <int>,
        "total_links": <int>,
        "total_clusters": <int>
      },
      "nodes": [
        {
          "id":           <str>,   // sample name
          "num_sequences":<int>,   // number of haplotype sequences
          "cluster":      <int|null>
        }, ...
      ],
      "links": [
        {
          "source":           <str>,
          "target":           <str>,
          "distance":         <float>,  // minimum_distance (alias for D3)
          "num_seqs_source":  <int>,
          "num_seqs_target":  <int>
        }, ...
      ]
    }
    """
    nodes = []
    for sid in sequences_dict:
        nodes.append({
            "id": sid,
            "num_sequences": len(sequences_dict[sid]),
            "cluster": sample_to_cluster.get(sid, None)
        })

    link_records = []
    for _, row in links_df.iterrows():
        link_records.append({
            "source": row['source'],
            "target": row['target'],
            "distance": row['minimum_distance'],       # field name D3 uses
            "num_seqs_source": int(row['num_seqs_source']),
            "num_seqs_target": int(row['num_seqs_target'])
        })

    num_clusters = len(set(sample_to_cluster.values())) if sample_to_cluster else 0

    network = {
        "metadata": {
            "threshold": threshold,
            "total_samples": len(nodes),
            "total_links": len(link_records),
            "total_clusters": num_clusters
        },
        "nodes": nodes,
        "links": link_records
    }

    return network


def create_sample_fasta():
    """Create a minimal sample FASTA file (no frequency in headers)."""
    sample_fasta = """\
>Sample_001_seq_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample_001_seq_2
ATCGATCGATCGAACGATCGATCGATCGATCGATCGATCG
>Sample_002_seq_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample_002_seq_2
ATCGATCGATCGAACGATCGATCGATCGATCGATCGATCG
>Sample_003_seq_1
GGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample_003_seq_2
GGCGATCGATCGTTCGATCGATCGATCGATCGATCGATCG
"""
    with open('sample_sequences.fasta', 'w') as f:
        f.write(sample_fasta)
    print("✅ Created sample file: sample_sequences.fasta")


def main():
    parser = argparse.ArgumentParser(
        description='Generate GHOST transmission network files from FASTA sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
FASTA file format:
  >SampleID_seqN
  SEQUENCE

  Example:
  >Sample_001_seq_1
  ATCGATCGATCG...

Examples:
  python real_data_processor.py --input sequences.fasta
  python real_data_processor.py --input sequences.fasta --threshold 0.05 --output ./results
  python real_data_processor.py --create-sample
        """
    )

    parser.add_argument('--input', '-i', type=str, help='Input FASTA file')
    parser.add_argument('--threshold', '-t', type=float, default=0.037,
                        help='Transmission distance threshold (default: 0.037 for HCV).If HAV, should be 0.')
    parser.add_argument('--output', '-o', type=str, default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--create-sample', action='store_true',
                        help='Create a sample FASTA file and exit')

    args = parser.parse_args()

    if args.create_sample:
        create_sample_fasta()
        return

    if not args.input:
        parser.print_help()
        return

    print("=" * 70)
    print("Real Data Processor")
    print("=" * 70)
    print(f"Input file : {args.input}")
    print(f"Threshold  : {args.threshold}")
    print(f"Output dir : {args.output}")
    print("=" * 70)

    # Parse FASTA
    print("\nParsing FASTA file...")
    sequences_dict = parse_fasta(args.input)
    if not sequences_dict:
        print("❌ Failed to parse FASTA file.")
        return
    print(f"✅ Parsed {len(sequences_dict)} samples, "
          f"{sum(len(v) for v in sequences_dict.values())} sequences total")

    # Calculate distances and build links
    print("Calculating pairwise Hamming distances...")
    links_df = calculate_pairwise_distances(sequences_dict, args.threshold)

    # Identify transmission clusters
    print("Identifying transmission clusters...")
    all_sample_ids = list(sequences_dict.keys())
    sample_to_cluster = identify_clusters_from_links(links_df, all_sample_ids)

    # Build cluster summary DataFrame
    cluster_rows = []
    for sid in all_sample_ids:
        cid = sample_to_cluster.get(sid, None)
        cluster_rows.append({
            'sample': sid,
            'cluster': cid if cid is not None else '',
            'num_sequences': len(sequences_dict[sid]),
        })
    cluster_df = pd.DataFrame(cluster_rows)

    # Annotate links with cluster ID
    if not links_df.empty:
        links_df['cluster'] = links_df['source'].map(
            lambda s: sample_to_cluster.get(s, '')
        )

    # Build network JSON
    print("Exporting network JSON...")
    network = export_network_json(sequences_dict, links_df, sample_to_cluster, args.threshold)

    # Save all outputs
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    links_file   = output_path / 'threshold_links.csv'
    cluster_file = output_path / 'cluster_summary.csv'
    json_file    = output_path / 'network.json'

    links_df.to_csv(links_file, index=False)
    cluster_df.to_csv(cluster_file, index=False)
    with open(json_file, 'w') as f:
        json.dump(network, f, indent=2)

    print(f"\n✅ Saved: {links_file}")
    print(f"✅ Saved: {cluster_file}")
    print(f"✅ Saved: {json_file}")

    # Summary stats
    linked_samples = (set(links_df['source']).union(set(links_df['target']))
                      if not links_df.empty else set())
    num_clusters = len(set(sample_to_cluster.values()))
    print(f"\nSummary:")
    print(f"  Samples              : {len(sequences_dict)}")
    print(f"  Transmission links   : {len(links_df)}")
    print(f"  Linked samples       : {len(linked_samples)}")
    print(f"  Unlinked samples     : {len(sequences_dict) - len(linked_samples)}")
    print(f"  Transmission clusters: {num_clusters}")
    print(f"\n💡 Open transmission-network.html and upload network.json to visualize.")


if __name__ == '__main__':
    main()
