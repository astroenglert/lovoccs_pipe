#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J final_visit_summary_process_band_cluster_name
#SBATCH -o slurm_outputs/final_visit_summary_process_band_cluster_name-%j.out
#SBATCH -e slurm_outputs/final_visit_summary_process_band_cluster_name-%j.err

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

# build the pipeline
pipetask build -p DRP-LoVoCCS.yaml#step2d \
    -C finalizeCharacterization:${CLUSTER_DIR}/configs/final_visit_summary_config_${BAND}.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step2d_${BAND}.yaml

# step2d final visit summary
# needed to add the calexp collection to include a few extra datasets... fortunately the pipeline is smart enough to only process visits which it has all the required datasets for
bps submit -b repo/repo \
    -i DECam/processing/quality_detectors_${BAND},DECam/processing/jointcal_${BAND},DECam/processing/visit_summary_${BAND},skymaps \
    -o DECam/processing/final_visit_summary_${BAND} \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step2d_${BAND}.yaml \
    -d "instrument='DECam' AND band='${BAND}' AND skymap='${CLN}_skymap'" \
    ${CLUSTER_DIR}/configs/bps_config.yaml

# technically we can do step2e as well... but its not required for deepCoaad
# if we ever want to do some difference-imaging for time-domain obsv, we'll want to do 2e


