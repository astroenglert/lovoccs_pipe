#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -J visit_summary_process_band_cluster_name
#SBATCH -o slurm_outputs/visit_summary_process_band_cluster_name-%j.out
#SBATCH -e slurm_outputs/visit_summary_process_band_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PY_SCRIPTS="py_scripts"
BAND=process_band

# navigate to .../cluster_name
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

# prevent implicit multithreading (otherwise tasks compete for resources)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# step2a visit summary
bps submit -b repo/repo \
    -i DECam/processing/quality_detectors_${BAND},skymaps \
    -o DECam/processing/visit_summary_${BAND} \
    -p DRP-LoVoCCS.yaml#step2a \
    -d "instrument='DECam' AND band='${BAND}'" \
    ${CLUSTER_DIR}/configs/bps_config.yaml


