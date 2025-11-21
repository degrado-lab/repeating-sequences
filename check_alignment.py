import os

def read_fasta(filepath):
    sequences = []
    with open(filepath, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                    seq = ''
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences

def write_fasta(sequences, output_path):
    with open(output_path, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f">rotated_seq_{i+1}\n{seq}\n")

def count_hydrophobics_per_position(sequences):
    hydrophobics = set("AVILMFYW")
    mer_length = len(sequences[0])
    counts = [0] * mer_length

    for seq in sequences:
        for i, aa in enumerate(seq):
            if aa in hydrophobics:
                counts[i] += 1

    return counts

def find_best_i_i2_pair(hydrophobic_counts):
    mer_length = len(hydrophobic_counts)
    best_i = None
    best_sum = -1

    for i in range(mer_length):
        i2 = (i + 2) % mer_length
        count_sum = hydrophobic_counts[i] + hydrophobic_counts[i2]
        if count_sum > best_sum:
            best_sum = count_sum
            best_i = i

    return best_i, (best_i + 2) % mer_length, best_sum

def rotate_sequence(seq, current_i, desired_i=0):
    mer_length = len(seq)
    shift = (current_i - desired_i) % mer_length
    return seq[shift:] + seq[:shift]

def run(config):
    output_dir = config["output_dir"]
    mer_length = config.get("repeat_unit_length", 6)

    clusters = [d for d in os.listdir(output_dir)
                if os.path.isdir(os.path.join(output_dir, d)) and d.startswith("cluster_")]

    for cluster in clusters:
        cluster_path = os.path.join(output_dir, cluster, "consensus_sequences")
        input_fasta = os.path.join(cluster_path, f"{cluster}_aligned_sequences.fasta")
        output_fasta = os.path.join(cluster_path, f"{cluster}_aligned_sequences_rotated.fasta")

        if not os.path.exists(input_fasta):
            print(f"⚠️ FASTA file not found: {input_fasta}")
            continue

        sequences = read_fasta(input_fasta)
        if not sequences:
            print(f"⚠️ No sequences found in {input_fasta}")
            continue

        hydrophobic_counts = count_hydrophobics_per_position(sequences)
        print(f"\n🧬 Hydrophobic counts per position in {cluster}:")
        for i, count in enumerate(hydrophobic_counts, start=1):
            print(f"  Position {i}: {count} hydrophobic residues")

        i, i2, total = find_best_i_i2_pair(hydrophobic_counts)
        print(f"🏆 Best i/i+2 hydrophobic pair in {cluster}:")
        print(f"  i = Position {i+1}, i+2 = Position {i2+1}, Total hydrophobic count = {total}")

        print(f"\n🔁 Preview rotation for cluster {cluster}:")
        print(f"  First sequence before: {sequences[0]}")
        print(f"  Rotating position {i+1} ({sequences[0][i]}) to front")

        rotated_sequences = [rotate_sequence(seq, i) for seq in sequences]

        print(f"  First sequence after : {rotated_sequences[0]}")

        write_fasta(rotated_sequences, output_fasta)
        print(f"✅ Rotated sequences written to {output_fasta}")

