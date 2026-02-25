#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J subtract_stars_cluster_name
#SBATCH -o slurm_outputs/subtract_stars_cluster_name-%j.out
#SBATCH -e slurm_outputs/subtract_stars_cluster_name-%j.err

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

# build the pipeline with our custom-configs
#TODO implement into central pipeline
pipetask build -p test.yaml#test \
    -C processBrightStars:${CLUSTER_DIR}/configs/subtract_stars_process_config.py \
    -C measureExtendedPsf:${CLUSTER_DIR}/configs/subtract_stars_measure_config.py \
    -C subtractBrightStars:${CLUSTER_DIR}/configs/subtract_stars_subtract_config.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/test_stars.yaml

# skycorr
#TODO pass the good selected visits, raws, and skyframes
bps submit -b repo/repo \
    -i DECam/processing/quality_detectors_g,DECam/processing/quality_detectors_r,DECam/processing/quality_detectors_i,DECam/processing/quality_detectors_z,DECam/calib/skyframes,DECam/processing/skycorr,refcats \
    -o DECam/processing/subtract_stars \
    -p ${CLUSTER_DIR}/pipeline_yamls/test_stars.yaml \
    -d "instrument='DECam' AND band IN ( 'g','r','i','z' )" \
    ${CLUSTER_DIR}/configs/bps_config.yaml

