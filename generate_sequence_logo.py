import os
import re
import logomaker
import pandas as pd
import matplotlib.pyplot as plt

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

def get_custom_color_scheme():
    return {
        'R': '#FF5500', 'K': '#FF5500', 'D': '#FF5500', 'E': '#FF5500',  # acidic/basic
        'Q': '#FF8844', 'N': '#FF8844', 'H': '#FF8844',                  # polar
        'S': '#FFBB88', 'T': '#FFBB88', 'C': '#FFBB88',                  # small polar
        'A': '#88BBFF',                                                  # small hydrophobic
        'Y': '#4488FF', 'W': '#4488FF', 'F': '#4488FF',                  # aromatics
        'I': '#0044FF', 'V': '#0044FF', 'L': '#0044FF', 'M': '#0044FF',  # aliphatic
        'G': '#000000',                                                  # glycine
        'P': '#FFD700',                                                  # proline
    }

def rotate_list(lst, offset):
    return lst[offset:] + lst[:offset]

def classify_apble(resname, phi, psi):
    if resname == "PRO":
        return "P"
    try:
        phi, psi = float(phi), float(psi)
    except:
        return "L"

    # α-helix region
    if -90 <= phi <= -30 and -80 <= psi <= -10:
        return "A"
    # β-sheet region
    elif (-180 <= phi <= -90 and 90 <= psi <= 180):
        return "B"
    # left-handed / extended region
    elif 30 <= phi <= 90 and 0 <= psi <= 90:
        return "E"
    else:
        return "L"

def run(config):
    output_root = config["output_dir"]
    color_scheme = get_custom_color_scheme()

    clusters = sorted([
        d for d in os.listdir(output_root)
        if d.startswith("cluster_") and os.path.isdir(os.path.join(output_root, d))
    ])

    for cluster in clusters:
        cluster_dir = os.path.join(output_root, cluster)
        fasta_path = os.path.join(cluster_dir, "consensus_sequences", f"{cluster}_aligned_sequences_rotated.fasta")
        phi_psi_path = os.path.join(cluster_dir, "phi_psi", f"{cluster}_phi_psi.txt")
        log_path = os.path.join(cluster_dir, "check_alignment_helix.log")
        logo_dir = os.path.join(cluster_dir, "sequence_logos")
        os.makedirs(logo_dir, exist_ok=True)

        print(f"🔍 Checking files for {cluster}")

        if not os.path.isfile(fasta_path) or not os.path.isfile(phi_psi_path):
            print(f"⚠️ Missing input files for {cluster}")
            continue

        sequences = read_fasta(fasta_path)
        if not sequences:
            print(f"⚠️ No sequences found in {fasta_path}")
            continue

        # Read rotation shift
        best_i = 0
        if os.path.exists(log_path):
            with open(log_path) as f:
                for line in f:
                    if "rotation_shift" in line:
                        best_i = int(line.strip().split(":")[1])
                        break

        # Parse phi/psi → APBLE
        apble_labels = []
        with open(phi_psi_path) as f:
            for line in f:
                match = re.match(
                    r"(\w+)-\d+:\s*\(\s*([-+]?\d*\.?\d+),\s*([-+]?\d*\.?\d+)\s*\)",
                    line.strip()
                )
                if match:
                    resname, phi, psi = match.groups()
                    label = classify_apble(resname, phi, psi)
                    apble_labels.append(label)
                else:
                    apble_labels.append("L")

        apble_labels = rotate_list(apble_labels, best_i)

        # Make logo
        df = logomaker.alignment_to_matrix(sequences, to_type='information')
        plt.figure(figsize=(4, 12))
        logo = logomaker.Logo(df, color_scheme=color_scheme)
        logo.ax.set_xticks(range(len(df)))
        logo.ax.set_xticklabels([str(i+1) for i in range(len(df))])
        logo.ax.set_title(f"{cluster}", fontsize=14, pad=30)
        logo.ax.set_ylabel("bits", fontsize=12)
        logo.ax.set_xlabel("Register Position", fontsize=12)
        logo.ax.spines['right'].set_visible(False)
        logo.ax.spines['top'].set_visible(False)

        # Annotate APBLE above sequence logo — using axes transform
        for i, label in enumerate(apble_labels[:len(df)]):
            logo.ax.text(
                i, 1.10, label,  # 1.10 = 10% above top of axis
                fontsize=12, fontweight='bold',
                ha='center', va='bottom',
                color='black',
                transform=logo.ax.get_xaxis_transform()
            )

        plt.tight_layout()
        output_file = os.path.join(logo_dir, f"{cluster}_sequence_logo.pdf")
        plt.savefig(output_file, dpi=300)
        plt.close()
        print(f"🧬 Saved: {output_file}")
