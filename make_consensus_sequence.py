import os
from collections import Counter

def read_fasta_sequences(filepath):
    sequences = []
    with open(filepath, 'r') as f:
        current_seq = ""
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line.strip()
        if current_seq:
            sequences.append(current_seq)
    return sequences

def compute_consensus(sequences):
    if not sequences:
        return ""

    length = len(sequences[0])
    consensus = []

    for i in range(length):
        column = [seq[i] for seq in sequences if len(seq) > i]
        most_common = Counter(column).most_common(1)[0][0]
        consensus.append(most_common)

    return ''.join(consensus)

def run(config):
    output_root = config["output_dir"]

    for cluster in os.listdir(output_root):
        if not cluster.startswith("cluster_"):
            continue

        cluster_dir = os.path.join(output_root, cluster)
        consensus_dir = os.path.join(cluster_dir, "consensus_sequences")
        os.makedirs(consensus_dir, exist_ok=True)

        aligned_fasta = os.path.join(consensus_dir, f"{cluster}_aligned_sequences_rotated.fasta")
        if not os.path.isfile(aligned_fasta):
            print(f"⚠️ Missing aligned sequences for {cluster}")
            continue

        sequences = read_fasta_sequences(aligned_fasta)
        if not sequences:
            print(f"⚠️ No sequences found in {aligned_fasta}")
            continue

        consensus_seq = compute_consensus(sequences)
        out_path = os.path.join(consensus_dir, f"{cluster}_consensus_sequence.fasta")

        with open(out_path, 'w') as f:
            f.write(f">{cluster}_consensus\n{consensus_seq}\n")

        print(f"✅ {cluster}: Consensus sequence saved to {out_path}")

