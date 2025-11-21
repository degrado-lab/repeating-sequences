import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def run(config):
    project_root = config["project_root"]
    output_dir = config["output_dir"]

    # Derive mer size from the project_root folder name (e.g., "5mer")
    mer_size = os.path.basename(project_root.rstrip("/"))

    csv_path = os.path.join(project_root, f"cluster_counts_{mer_size}.csv")
    output_pdf = os.path.join(output_dir, f"{mer_size}_top_clusters_histogram.pdf")
  
    print(f"🔍 Reading: {csv_path}")
    print(f"📊 Will save to: {output_pdf}")

    if not os.path.isfile(csv_path):
        print(f"⚠️ Missing: {csv_path}")
        return

    try:
        df = pd.read_csv(csv_path)
        print(f"📄 Columns: {df.columns.tolist()}")
    except Exception as e:
        print(f"❌ Failed to read {csv_path}: {e}")
        return

    if 'cluster_label' not in df.columns or 'count' not in df.columns:
        print(f"⚠️ Required columns not found in {csv_path}")
        return

    # Sort and select top 10 (or fewer)
    df_sorted = df.sort_values(by="count", ascending=False).head(10)
    num_clusters = len(df_sorted)

    if num_clusters == 0:
        print("⚠️ No clusters to plot.")
        return

    # --- Prepare data ---
    x = np.arange(num_clusters)
    labels = df_sorted["cluster_label"].astype(str)
    counts = df_sorted["count"]

    # --- Plot setup ---
    fig, ax = plt.subplots(figsize=(20, 5))
    ax.bar(x, counts, color="#1f77b4", width=0.9)

    # --- Add a slight margin so bars aren't flush ---
    pad = 0.6  # tweak between 0.3–0.6 if needed
    ax.set_xlim(-pad, num_clusters - 1 + pad)

    # --- Axis labels & ticks ---
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0)
    ax.set_ylabel("Count", fontsize=18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    # --- Clean up spines ---
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.spines['left'].set_linewidth(1.2)

    # --- Y-axis limit (dynamic) ---
    ax.set_ylim(0, max(counts) * 1.1)

    # --- Layout and export ---
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    fig.savefig(output_pdf, bbox_inches="tight", transparent=True)
    plt.close(fig)
    plt.show()

    print(f"✅ Saved histogram to: {output_pdf}")

