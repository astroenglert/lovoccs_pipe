#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -J export_data_cluster_name
#SBATCH -o slurm_outputs/export_data_cluster_name-%j.out
#SBATCH -e slurm_outputs/export_data_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PY_SCRIPTS="py_scripts"

# navigate to .../A85
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

python -m python_scripts.export.export_data "${CLN}"


