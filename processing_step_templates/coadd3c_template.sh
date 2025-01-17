#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J coadd_3c_cluster_name
#SBATCH -o slurm_outputs/coadd_3c_cluster_name-%j.out
#SBATCH -e slurm_outputs/coadd_3c_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER

# navigate to .../cluster_name
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# prevent implicit multithreading (otherwise tasks compete for resources)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# build the pipeline
pipetask build -p DRP-LoVoCCS.yaml#step3c \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step3c.yaml

# step3a assembleCoadd, detect, deblend
bps submit -b repo/repo \
    -i DECam/processing/coadd_3a,DECam/processing/coadd_3b_u,DECam/processing/coadd_3b_g,DECam/processing/coadd_3b_r,DECam/processing/coadd_3b_i,DECam/processing/coadd_3b_z,skymaps \
    -o DECam/processing/coadd_3c \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step3c.yaml \
    -d "instrument='DECam' AND skymap='${CLN}_skymap'" \
    ${CLUSTER_DIR}/configs/bps_config.yaml


