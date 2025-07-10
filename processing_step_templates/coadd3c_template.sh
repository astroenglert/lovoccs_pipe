#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J coadd_3d_cluster_name
#SBATCH -o slurm_outputs/coadd_3d_cluster_name-%j.out
#SBATCH -e slurm_outputs/coadd_3d_cluster_name-%j.err

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

# prevent implicit multithreading (otherwise tasks compete for resources)
export OMP_NUM_THREADS=1

# build the pipeline with our custom-configs
pipetask build -p DRP-LoVoCCS.yaml#step3c \
    -C skyCorr:${CLUSTER_DIR}/configs/coadd_3c_skycorr_config.py \
    -c makeWarp:doApplySkyCorr=True \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step3c.yaml

# skycorr
bps submit -b repo/repo \
    -i DECam/processing/quality_detectors_g,DECam/processing/quality_detectors_r,DECam/processing/quality_detectors_i,DECam/processing/quality_detectors_z,DECam/processing/final_visit_summary_g,DECam/processing/final_visit_summary_r,DECam/processing/final_visit_summary_i,DECam/processing/final_visit_summary_z,DECam/calib/certified,DECam/calib/unbounded,DECam/raw/all \
    -o DECam/processing/skycorr \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step3c.yaml \
    -d "instrument='DECam' AND band IN ( 'g','r','i','z' )" \
    ${CLUSTER_DIR}/configs/bps_config_coadd3a.yaml

