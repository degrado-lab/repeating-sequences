import os
import re

def parse_tmalign_output(tmalign_output, mer_length):
    sequences = []
    current_top = None
    current_align = None
    current_bot = None
    match_count = 0
    target_seq = None

    for line in tmalign_output.splitlines():
        line = line.strip()

        if re.match(r'^[A-Z\-]+$', line):
            if current_top is None:
                current_top = line
            elif current_bot is None:
                current_bot = line
        elif re.match(r'^[.: ]+$', line):
            current_align = line

        if current_top and current_align and current_bot:
            print(f"\n🧩 Found alignment block:")
            print(f"TOP: {current_top}")
            print(f"ALI: {current_align}")
            print(f"BOT: {current_bot}")

            if target_seq is None:
                target_seq = ''.join([c for c in current_bot if c != '-'])[:mer_length]
                print(f"🎯 Target sequence: {target_seq}")

            match_seqs = []
            for i in range(len(current_bot) - mer_length + 1):
                window = current_bot[i:i + mer_length]
                if window == target_seq:
                    matched = current_top[i:i + mer_length]
                    match_seqs.append(matched)

            print(f"🔍 Found {len(match_seqs)} matches in BOT line.")

            for idx, seq in enumerate(match_seqs[1:], start=2):
                print(f"  Checking match #{idx}: {seq}")
                if '-' not in seq:
                    sequences.append(seq)
                    match_count += 1
                    print(f"✅ Accepted: {seq}")
                    break
                else:
                    print(f"⛔ Rejected (gap): {seq}")

            current_top = current_align = current_bot = None

    print(f"\n🔎 Total accepted sequences: {len(sequences)}")
    return sequences

def write_fasta(sequences, output_fasta):
    with open(output_fasta, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq_{i+1}\n{seq}\n")

def count_hydrophobics_per_position(sequences):
    hydrophobics = set("AVILMFYW")  # Classical hydrophobic residues
    mer_length = len(sequences[0])
    counts = [0] * mer_length

    for seq in sequences:
        for i, aa in enumerate(seq):
            if aa in hydrophobics:
                counts[i] += 1

    return counts

def run(config):
    input_dir = config["input_dir"]
    output_dir = config["output_dir"]
    mer_length = config.get("repeat_unit_length", 6)

    clusters = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d)) and d.startswith("cluster_")]

    for cluster in clusters:
        cluster_output_dir = os.path.join(output_dir, cluster)
        consensus_dir = os.path.join(cluster_output_dir, "consensus_sequences")
        os.makedirs(consensus_dir, exist_ok=True)

        input_file = os.path.join(cluster_output_dir, "clusters_to_members.txt")
        output_fasta = os.path.join(consensus_dir, f"{cluster}_aligned_sequences.fasta")

        print(f"\n🔹 Parsing TM-align output for aligned sequences in {cluster}...")

        if not os.path.exists(input_file):
            print(f"⚠️ TM-align file not found: {input_file}")
            continue

        with open(input_file, 'r') as file:
            tmalign_text = file.read()

        sequences = parse_tmalign_output(tmalign_text, mer_length)

        if sequences:
            write_fasta(sequences, output_fasta)
            print(f"✅ Written {len(sequences)} sequences for {cluster} to {output_fasta}")

            # 🧪 Count hydrophobics per position
            hydrophobic_counts = count_hydrophobics_per_position(sequences)
            print(f"🧬 Hydrophobic counts per position in {cluster}:")
            for i, count in enumerate(hydrophobic_counts, start=1):
                print(f"  Position {i}: {count} hydrophobic residues")

        else:
            print(f"⚠️ No suitable aligned sequences found for {cluster}")

