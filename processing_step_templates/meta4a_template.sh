#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J meta_4a_cluster_name
#SBATCH -o slurm_outputs/meta_4a_cluster_name-%j.out
#SBATCH -e slurm_outputs/meta_4a_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PY_SCRIPTS="py_scripts"

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
pipetask build -p DRP-LoVoCCS.yaml#step4a \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step4a.yaml

# only run metadetect on the central patches where data is deepest
# a little tighter than 33-88 since those corners are mostly empty data
PATCHES="(40..43, 51..56, 63..68, 75..80, 87..92, 100..103)"

# step3a assembleCoadd, detect, deblend
bps submit -b repo/repo \
    -i DECam/processing/coadd_3a,refcats,skymaps \
    -o DECam/processing/meta_4a \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step4a.yaml \
    -d "instrument='DECam' AND skymap='${CLN}_skymap' AND ( patch IN ${PATCHES} )" \
    ${CLUSTER_DIR}/configs/bps_config.yaml


