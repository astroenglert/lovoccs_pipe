#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J skycorr_cluster_name
#SBATCH -o slurm_outputs/skycorr_cluster_name-%j.out
#SBATCH -e slurm_outputs/skycorr_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="/gpfs/data/idellant/englert_newPipelineDev/A85" # UPDATE LATER
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

# first, run skycorr to create sky-corrected backgrounds
# build the pipeline with our custom-configs
#TODO centralize our custom version of DRP-Merian
pipetask build -p DRP-Merian_copy.yaml#stepSky \
    -C skyCorr:${CLUSTER_DIR}/configs/skycorr_config.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_skycorr.yaml

# skycorr
bps submit -b repo/repo \
    -i DECam/processing/quality_detectors_g,DECam/processing/quality_detectors_r,DECam/processing/quality_detectors_i,DECam/processing/quality_detectors_z,DECam/calib/skyframes,DECam/calib/unbounded,DECam/raw/all \
    -o DECam/processing/skycorr \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_skycorr.yaml \
    -d "instrument='DECam' AND band IN ( 'g','r','i','z' )" \
    ${CLUSTER_DIR}/configs/bps_config.yaml

# now run coaddition with those, we simply re-run coadd_3a with updated configs for this part
#TODO test catalog-level products derived from skycorr coadded frames, see notes in detection_config_template

# build the pipeline
pipetask build -p DRP-Merian_copy.yaml#step3a \
    -C assembleCoadd:configs/skycorr_assemblecoadd_config.py \
    -C makeWarp:configs/skycorr_makewarp_config.py \
    -C detection:configs/skycorr_detection_config.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step3a_skycorr.yaml

# step3a assembleCoadd, detect, deblend
bps submit -b repo/repo \
    -i DECam/processing/quality_detectors_g,DECam/processing/quality_detectors_r,DECam/processing/quality_detectors_i,DECam/processing/quality_detectors_z,\
DECam/processing/final_visit_summary_g,DECam/processing/final_visit_summary_r,DECam/processing/final_visit_summary_i,DECam/processing/final_visit_summary_z,\
skymaps,DECam/processing/skycorr \
    -o DECam/processing/coadd_skycorr \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step3a_skycorr.yaml \
    -d "instrument='DECam' AND skymap='${CLN}_skymap' AND band IN ( 'g','r','i','z' )" \
    ${CLUSTER_DIR}/configs/bps_config.yaml

