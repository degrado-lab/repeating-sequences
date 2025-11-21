import os
import subprocess
import re

def get_clusters(config):
    cluster_dir = config["all_clusters_dir"]
    return sorted([
        os.path.splitext(f)[0]
        for f in os.listdir(cluster_dir)
        if f.endswith(".pdb")
    ])

def parse_rmsd(output):
    """Extract RMSD value from TM-align output."""
    match = re.search(r"RMSD=(\s*\d+\.\d+)", output)
    if match:
        return float(match.group(1).strip())
    return None

def run(config):
    all_clusters_dir = config["all_clusters_dir"]
    output_dir = config["output_dir"]
    tmalign_path = config["tmalign_path"]
    rmsd_cutoff = config.get("clusterwise_rmsd_cutoff", 1.8)
    mer_len = config.get("repeat_unit_length", 5)
    trim_len = mer_len * 6  # number of residues to include from N-term

    clusters = get_clusters(config)

    os.makedirs(output_dir, exist_ok=True)

    for i, cluster_a in enumerate(clusters):
        cluster_a_path = os.path.join(all_clusters_dir, f"{cluster_a}.pdb")
        cluster_output_dir = os.path.join(output_dir, cluster_a)
        os.makedirs(cluster_output_dir, exist_ok=True)

        all_output_file = os.path.join(cluster_output_dir, f"{cluster_a}_vs_clusters.txt")
        low_rmsd_file = os.path.join(cluster_output_dir, f"{cluster_a}_below_{rmsd_cutoff}_rmsd.txt")

        with open(all_output_file, "w") as all_f, open(low_rmsd_file, "w") as low_f:
            for cluster_b in clusters:
                if cluster_a == cluster_b:
                    continue

                cluster_b_path = os.path.join(all_clusters_dir, f"{cluster_b}.pdb")

                cmd = [
                    tmalign_path,
                    cluster_a_path,
                    cluster_b_path,
                    "-ter",
                    str(trim_len)
                ]

                try:
                    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
                    output = result.stdout

                    all_f.write(f">>> {cluster_a} vs {cluster_b}\n")
                    all_f.write(output + "\n\n")

                    rmsd = parse_rmsd(output)
                    if rmsd is not None and rmsd < rmsd_cutoff:
                        low_f.write(f">>> {cluster_a} vs {cluster_b} — RMSD = {rmsd:.2f}\n")
                        low_f.write(output + "\n\n")

                except subprocess.CalledProcessError:
                    print(f"❌ TM-align failed: {cluster_a} vs {cluster_b}")
                    continue

        print(f"✅ Finished comparing {cluster_a} to all clusters → {all_output_file}")

