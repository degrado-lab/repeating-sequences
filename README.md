**Repeating Sequence Cluster Analysis & TM-align Search Pipelines**
This repository contains two major workflows:
A full cluster-analysis pipeline for repeated protein units (repeat units, solenoids, etc.).
A TM-align–based database-search pipeline for screening a structural database (e.g., PDB60) against a target repeat structure.
Both pipelines are modular, configurable, and designed to be run end-to-end.

**1. Cluster Analysis Pipeline**
The cluster analysis workflow organizes structural clusters, aligns structural repeats, filters based on RMSD/TM-align scores, generates consensus sequences, and creates diagnostic plots and logos.

**Repository Structure**
repo/
│
├── run_pipeline.py
├── config.yaml
├── helper_scripts/
│   ├── organize_clusters.py
│   ├── run_tmalign_predictions.py
│   ├── cluster_ct_histograms.py
│   ├── clusterwise_tmalign_filter.py
│   ├── align_sequences.py
│   ├── check_alignment.py
│   ├── make_consensus_sequence.py
│   ├── generate_sequence_logo.py
│   └── cp_phi_psi.py
│
└── tmalign_database_search/
    ├── tmalign.sh
    └── high_tmscore.py

**Config File (config.yaml)**
You control nearly all pipeline behavior using config.yaml.
Example:

cluster_size: 10
repeat_unit_length: 10
project_root: "/path/to/your/repeat/folder"

tmalign_path: "path/to/your/TMalign"

clusterwise_rmsd_cutoff: 1.2
use_gpu: true

input_dir: "{project_root}/all_pdbs"
output_dir: "{project_root}/outputs"
all_clusters_dir: "{project_root}/all_clusters"
all_pdbs_dir: "{project_root}/all_pdbs"

cluster_csv_file: "{project_root}/cluster_counts_10mer.csv"
clusters_to_members_file: "{output_dir}/clusters_to_members.txt"

Any string values containing {project_root} (or other config keys) are auto-expanded by the pipeline.

**Running the Full Pipeline**
Simply execute:

python run_pipeline.py

The script automatically:

loads the config
resolves paths
detects clusters in all_clusters_dir
runs all helper scripts in sequence

**Helper Script Descriptions**
Each helper script corresponds to a step in the pipeline.
Below are concise descriptions suitable for documentation or internal notes.

1. organize_clusters.py
Organizes cluster member PDB files into per-cluster directories based on the cluster‐counts CSV. This prepares the folder structure used in all subsequent analysis.

2. run_tmalign_predictions.py
Runs TM-align between each cluster centroid and its member structures, storing alignment logs and scores for downstream filtering.

3. cluster_ct_histograms.py
Generates histograms summarizing repeat cluster sizes and distributions. Useful for quality-control and dataset overview.

4. clusterwise_tmalign_filter.py
Filters TM-align results based on RMSD/TM-score thresholds (e.g., clusterwise_rmsd_cutoff), keeping only well-aligned cluster members.

5. align_sequences.py
Builds sequence alignments for repeat units based on the structural TM-align correspondences. Produces per-cluster FASTA alignments or alignment text files.

6. check_alignment.py
Performs QC on the generated alignments (length checks, coverage checks, gap sanity, etc.) and reports clusters requiring manual inspection.

7. make_consensus_sequence.py
Constructs consensus sequences for each cluster based on the filtered alignments.

8. generate_sequence_logo.py
Creates publication-ready sequence logos (e.g., via Logomaker or WebLogo) from the cluster alignments.

9. cp_phi_psi.py
Generates φ/ψ (phi/psi) torsion-angle plots for visualizing structural variation across repeat units.

**2. TM-align Database-Search Pipeline**

The second major workflow screens an input repeat structure against a structural database (e.g., PDB60) using TM-align.

This pipeline is located in:

tmalign_database_search/

**Requirements**
TM-align (compiled binary)
A structural database in .cif form (e.g., PDB60)
macOS or Linux with bash

**Configuring & Running tmalign.sh**
Edit the top of the script:

TMALIGN_BIN=/path/to/TMalign
TIMEOUT_CMD="timeout"   # macOS: use gtimeout from coreutils
TIME_LIMIT="5s"

TARGET_PDB="/path/to/target_structure.pdb"
CIF_DIR="/path/to/pdb60_database"
OUTPUT_DIR="/path/to/output_directory"

Then run: 

bash tmalign.sh

This script ecursively scans the CIF directory
Runs TM-align of the target PDB vs each CIF file
Logs timeouts and errors
Writes individual TM-align logs to OUTPUT_DIR
This produces hundreds or thousands of TM-align output files.

**Filtering High TM-score Hits: high_tmscore.py**

After running tmalign.sh, run:

python high_tmscore.py

This script parses all .txt TM-align logs in the current directory
Extracts TM-score (normalized by Chain_1 length)

Applies a user-defined cutoff:
THRESHOLD = 0.5

Outputs a CSV file:
high_tm_scores.csv

Containing:
filename	tm_score

**Example Workflow**
Build a target structure → e.g., cluster centroid PDB
Download PDB60 (or any CIF database)
Install TM-align
Edit paths in tmalign.sh
Run the database search
Run high_tmscore.py to find good structural matches
Interpret results & follow up with clustering pipeline if desired

**Citation**
If you use this pipeline in a publication, please cite TM-align:
Zhang & Skolnick (2005). TM-align: A protein structure alignment algorithm based on the TM-score.
If you'd like a formal citation block for this repo, I can generate one.

**Contact**
If you encounter issues, feel free to open a GitHub Issue or email rose.yang@ucsf.edu.
