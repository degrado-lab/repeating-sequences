import os
import yaml

from helper_scripts import organize_clusters
from helper_scripts import run_tmalign_predictions
from helper_scripts import cluster_ct_histograms
from helper_scripts import clusterwise_tmalign_filter

from helper_scripts import align_sequences
from helper_scripts import check_alignment

from helper_scripts import make_consensus_sequence

from helper_scripts import generate_sequence_logo

from helper_scripts import cp_phi_psi



# LOAD CONFIGURATION

def load_config():
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)

    if "project_root" in config:
        for k, v in config.items():
            if isinstance(v, str):
                config[k] = v.format(**config)
            elif isinstance(v, list):
                config[k] = [x.format(**config) if isinstance(x, str) else x for x in v]

    cluster_root = config["all_clusters_dir"]
    if os.path.isdir(cluster_root):
        clusters = []
        for d in os.listdir(cluster_root):
            if d.startswith("cluster_") and (
                d.endswith(".pdb") or os.path.isdir(os.path.join(cluster_root, d))
            ):
                clusters.append(os.path.splitext(d)[0])
        clusters.sort()
        config["clusters_to_analyze"] = clusters
    else:
        print(f"⚠️ Warning: Cluster directory {cluster_root} not found.")
        config["clusters_to_analyze"] = []

    return config


# MAIN PIPELINE

def main():
    config = load_config()
    print("🏁 Starting Full Cluster Analysis Pipeline!")

    # Organize cluster member pdb files into cluster directories based on cluster_count csv file 
    print("\n📦 Organizing cluster outputs...")
    organize_clusters.run(config)

    # Run tmalign_prediction between cluster centroid and cluster member pdbs
    print("\n🔍 Running TM-align predictions (cluster → members)...")
    run_tmalign_predictions.run(config)

    print("\n📊 Generating cluster count histograms...")
    cluster_ct_histograms.run(config)

    print("\n🧹 Filtering clusterwise TM-align results...")
    clusterwise_tmalign_filter.run(config)

    # Aligning sequences of specified repeat length by their TM-align alignments
    print("\n🧬 Aligning sequences...")
    align_sequences.run(config)

    # Checking alignments 
    print("\n🧪 Checking alignments...")
    check_alignment.run(config)

    # Bulding consensus sequence from sequence alignment text files 
    print("\n🧬 Building consensus sequences...")
    make_consensus_sequence.run(config)

    # Generating phi_psi plots 
    print("\n🎨 Generating φ/ψ angle plots...")
    cp_phi_psi.run(config)

    # Generating sequence logos from sequence alignments 
    print("\n🔠 Generating sequence logos...")
    generate_sequence_logo.run(config)

    print("\n🎉 Pipeline complete!")


if __name__ == "__main__":
    main()

