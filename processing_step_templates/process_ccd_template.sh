#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J process_ccd_process_band_cluster_name
#SBATCH -o slurm_outputs/process_ccd_process_band_cluster_name-%j.out
#SBATCH -e slurm_outputs/process_ccd_process_band_cluster_name-%j.err

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

# generate xtalksources 
# before generating these, check if they've been created (process_ccd can be troublesome and recreating the xtalksources w. every iteration can lead to issues

# dump the query into a variable to read
OUT="$(butler query-collections repo/repo DECam/processing/xtalksources_${BAND})"

# if the chained collection does not exist, then compute the crosstalk sources
# bash is annoying for this type of thing and the cln-interface for LSP isn't meant to be used this way...
# but the alternative is writing a dedicated python script to run all this
if ! [[ ${OUT} == *"DECam/processing/xtalksources_${BAND}"*"CHAINED"* ]]; then

bps submit -b repo/repo \
    -i DECam/raw/all,\
DECam/calib/curated,DECam/calib/certified,DECam/calib/curated/curated/19700101T000000Z,DECam/calib/curated/unbounded \
    -o DECam/processing/xtalksources_${BAND} \
    -p DRP-LoVoCCS.yaml#step0 \
    -d "instrument='DECam' AND band='${BAND}'" \
    ${CLUSTER_DIR}/configs/bps_config.yaml

fi

# build the pipeline with our custom-configs
pipetask build -p DRP-LoVoCCS.yaml#step1 \
    -C calibrate:${CLUSTER_DIR}/configs/process_ccd_refcat_config_${BAND}.py \
    -C characterizeImage:${CLUSTER_DIR}/configs/process_ccd_psf_config_${BAND}.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step1_${BAND}.yaml

# register a dataset-type (bps won't do this automatically unlike pipetasks)
butler register-dataset-type repo/repo icSrc_schema SourceCatalog

# process_ccd
bps submit -b repo/repo \
    -i DECam/raw/all,\
DECam/calib/curated,DECam/calib/certified,DECam/calib/curated/curated/19700101T000000Z,DECam/calib/curated/unbounded,\
DECam/processing/xtalksources_${BAND},refcats \
    -o DECam/processing/calexp_${BAND} \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step1_${BAND}.yaml \
    -d "instrument='DECam' AND band='${BAND}'" \
    ${CLUSTER_DIR}/configs/bps_config.yaml

