import os
import subprocess

def count_residues(pdb_file):
    residues = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15].strip() == 'CA':
                resid = (line[21], line[22:26].strip())
                residues.add(resid)
    return len(residues)

def crop_pdb(input_pdb, output_pdb, max_residues):
    residues_seen = set()
    with open(input_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
        for line in fin:
            if line.startswith("ATOM") and line[13:15].strip() == 'CA':
                resid = (line[21], line[22:26].strip())
                if resid not in residues_seen:
                    residues_seen.add(resid)
                if len(residues_seen) > max_residues:
                    break
            fout.write(line)

def run(config):
    cluster_dir = config['all_clusters_dir']
    input_dir = config['input_dir']
    output_dir = config['output_dir']
    clusters_to_process = config['clusters_to_analyze']
    tmalign_path = config['tmalign_path']

    for cluster_name in clusters_to_process:
        cluster_output_dir = os.path.join(output_dir, cluster_name)
        os.makedirs(cluster_output_dir, exist_ok=True)

        cluster_center_pdb = os.path.join(cluster_dir, f"{cluster_name}.pdb")
        cluster_member_dir = os.path.join(input_dir, cluster_name)

        if not os.path.isfile(cluster_center_pdb):
            print(f"❌ Cluster center {cluster_center_pdb} not found!")
            continue
        if not os.path.isdir(cluster_member_dir):
            print(f"❌ Cluster member directory {cluster_member_dir} not found!")
            continue

        max_residues = count_residues(cluster_center_pdb)
        print(f"🔹 Processing {cluster_name} (Cluster center has {max_residues} residues)")

        output_file = os.path.join(cluster_output_dir, "clusters_to_members.txt")

        pdb_files = [f for f in os.listdir(cluster_member_dir) if f.endswith(".pdb")]

        with open(output_file, "w") as out:
            for pdb_file in pdb_files:
                member_pdb_path = os.path.join(cluster_member_dir, pdb_file)
                cropped_pdb_path = os.path.join(cluster_output_dir, f"cropped_{pdb_file}")

                crop_pdb(member_pdb_path, cropped_pdb_path, max_residues)

                try:
                    result = subprocess.run([tmalign_path, cropped_pdb_path, cluster_center_pdb],
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                            universal_newlines=True, check=True)
                    out.write(f"Comparing {pdb_file} to {cluster_name}.pdb\n")
                    out.write(result.stdout)
                    out.write("\n\n")
                except subprocess.CalledProcessError as e:
                    print(f"⚠️ TM-align failed for {pdb_file}: {e}")

                os.remove(cropped_pdb_path)

        print(f"✅ TM-align output saved: {output_file}")

