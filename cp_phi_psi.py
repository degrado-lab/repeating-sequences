import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import traceback

def load_apble_background():
    csv_path = "/Users/rose/Downloads/APBLE_Mode_36x36.csv"
    apble_df = pd.read_csv(csv_path, header=None)
    apble_colors = {
        'A': '#D8BFD8', 'P': '#90EE90', 'B': '#ADD8E6',
        'L': '#FFDAB9', 'E': '#FFFACD'
    }

    n = 36
    rgb_image = np.ones((n, n, 3))
    for i in range(n):
        for j in range(n):
            entry = apble_df.iat[i, j]
            letter = str(entry).strip() if pd.notna(entry) else ''
            if letter not in apble_colors:
                rgb_image[i, j, :] = (0.95, 0.95, 0.95)
            else:
                hex_color = apble_colors[letter].lstrip("#")
                rgb = tuple(int(hex_color[k:k+2], 16)/255 for k in (0, 2, 4))
                rgb_image[i, j, :] = rgb

    wrapped_rows = rgb_image[-2:, :, :]
    extended_rgb_image = np.vstack([rgb_image, wrapped_rows])
    return np.flipud(extended_rgb_image)

def generate_cyclic_permutations(seq):
    return [seq[i:] + seq[:i] for i in range(len(seq))]

def get_repeat_from_cluster_structure(full_sequence, repeat_len):
    return full_sequence[:repeat_len]

def match_rotated_repeat_in_fasta(permutations, fasta_path):
    with open(fasta_path, "r") as f:
        lines = [line.strip() for line in f if not line.startswith(">")]
    for perm in permutations:
        for line in lines:
            if perm in line:
                print(f"✅ Found rotated repeat '{perm}' in aligned sequences")
                return perm
    return None

def find_repeat_in_structure(all_residues, all_phi_psi, repeat_unit, aa3to1, occurrence=0):
    window_len = len(repeat_unit)
    matches = []

    for i in range(len(all_phi_psi) - window_len + 1):
        window_phi_psi = all_phi_psi[i:i + window_len]
        window_residues = all_residues[i:i + window_len]
        try:
            window_seq = "".join(aa3to1[res.get_resname()] for res in window_residues)
        except KeyError:
            continue
        if window_seq == repeat_unit:
            matches.append((window_residues, window_phi_psi))

    if len(matches) >= occurrence:
        return matches[occurrence - 1]  # 6th match
    else:
        print(f"⚠️ Only found {len(matches)} matches of {repeat_unit}, need {occurrence}")
        return None, None

def run(config):
    output_dir = config["output_dir"]
    repeat_len = config.get("repeat_unit_length", 5)
    apble_background = load_apble_background()

    aa3to1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
        'TRP': 'W', 'TYR': 'Y'
    }

    for cluster in config["clusters_to_analyze"]:
        pdb_path = os.path.join(config["all_clusters_dir"], f"{cluster}.pdb")
        if not os.path.isfile(pdb_path):
            print(f"❌ Missing PDB for {cluster}: {pdb_path}")
            continue

        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(cluster, pdb_path)
            ppb = PPBuilder()
            peptides = list(ppb.build_peptides(structure))
            if not peptides:
                print(f"⚠️ No peptide found in {cluster}")
                continue

            full_sequence = "".join(str(pep.get_sequence()) for pep in peptides)
            all_residues, all_phi_psi = [], []
            for pep in peptides:
                all_residues.extend(pep)
                all_phi_psi.extend(pep.get_phi_psi_list())

            cluster_repeat = get_repeat_from_cluster_structure(full_sequence, repeat_len)
            permutations = generate_cyclic_permutations(cluster_repeat)

            rotated_fasta_path = os.path.join(output_dir, cluster, "consensus_sequences", f"{cluster}_aligned_sequences_rotated.fasta")
            if not os.path.isfile(rotated_fasta_path):
                print(f"⚠️ Rotated FASTA not found: {rotated_fasta_path}")
                continue

            matched_repeat = match_rotated_repeat_in_fasta(permutations, rotated_fasta_path)
            if not matched_repeat:
                print(f"⚠️ No rotated match of {cluster_repeat} found in aligned sequences for {cluster}")
                continue

            subset_residues, subset_phi_psi = find_repeat_in_structure(
                all_residues, all_phi_psi, matched_repeat, aa3to1, occurrence=6
            )
            if subset_residues is None:
                print(f"⚠️ Could not find 6th occurrence of {matched_repeat} in structure")
                continue

            print(f"\n📐 φ/ψ angles for {cluster}:")
            print(f"🧬 Using 6th occurrence of rotated repeat: {matched_repeat}")
            phis, psis = [], []
            for res, (phi, psi) in zip(subset_residues, subset_phi_psi):
                resname = res.get_resname()
                resnum = res.get_id()[1]
                phi_deg = np.degrees(phi) if phi is not None else None
                psi_deg = np.degrees(psi) if psi is not None else None
                if psi_deg is not None and psi_deg < -160:
                    psi_deg += 360
                phi_str = f"{phi_deg:.1f}" if phi_deg is not None else "None"
                psi_str = f"{psi_deg:.1f}" if psi_deg is not None else "None"
                print(f"{resname}-{resnum}:  ( {phi_str}, {psi_str} )")
                if phi_deg is not None and psi_deg is not None:
                    phis.append(phi_deg)
                    psis.append(psi_deg)

            if len(phis) < 2:
                print(f"⚠️ Too few valid points to plot for {cluster}")
                continue

            y_lower, y_upper = (-180, 180)
            if any(psi > 180 for psi in psis):
                y_lower, y_upper = (-160, 200)

            out_dir = os.path.join(output_dir, cluster, "phi_psi")
            os.makedirs(out_dir, exist_ok=True)
            out_txt_path = os.path.join(out_dir, f"{cluster}_phi_psi.txt")
            with open(out_txt_path, "w") as out_txt:
                for res, (phi, psi) in zip(subset_residues, subset_phi_psi):
                    resname = res.get_resname()
                    resnum = res.get_id()[1]
                    phi_deg = np.degrees(phi) if phi is not None else None
                    psi_deg = np.degrees(psi) if psi is not None else None
                    if psi_deg is not None and psi_deg < -160:
                        psi_deg += 360
                    phi_str = f"{phi_deg:.1f}" if phi_deg is not None else "None"
                    psi_str = f"{psi_deg:.1f}" if psi_deg is not None else "None"
                    out_txt.write(f"{resname}-{resnum}:  ( {phi_str}, {psi_str} )\n")

            fig, ax = plt.subplots(figsize=(5, 5))
            ax.imshow(apble_background, extent=[-180, 180, -180, 200], origin='lower',
                      interpolation='nearest', aspect='auto', zorder=0)

            for i in range(len(phis)):
                x0, y0 = phis[i], psis[i]
                x1, y1 = phis[(i + 1) % len(phis)], psis[(i + 1) % len(phis)]
                ax.plot([x0, x1], [y0, y1], color='gray', linewidth=1.5, zorder=1)

            cmap = cm.get_cmap("Greys")
            colors = [cmap(i / (len(phis) - 1)) for i in range(len(phis))]
            for i, (phi, psi) in enumerate(zip(phis, psis)):
                ax.scatter(phi, psi, color=colors[i], s=40, edgecolor='black', zorder=2)

            ax.axhline(0, color="black", linewidth=1)
            ax.axvline(0, color="black", linewidth=1)
            ax.set_xlim(-180, 180)
            ax.set_ylim(y_lower, y_upper)
            ax.set_xlabel("Phi (°)")
            ax.set_ylabel("Psi (°)")
            ax.set_title(f"Ramachandran Plot: {cluster}")
            ax.set_aspect('equal')

            fig.subplots_adjust(left=0.18, right=0.95, bottom=0.12, top=0.90)
            out_path = os.path.join(out_dir, f"{cluster}_ramachandran.pdf")
            plt.savefig(out_path, dpi=300, pad_inches=0.05)
            plt.close()
            print(f"📈 Saved: {out_path}")

        except Exception as e:
            print(f"❌ Error for {cluster}: {e}")
            traceback.print_exc()

