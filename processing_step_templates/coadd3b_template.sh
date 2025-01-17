#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J coadd_3b_process_band_cluster_name
#SBATCH -o slurm_outputs/coadd_3b_process_band_cluster_name-%j.out
#SBATCH -e slurm_outputs/coadd_3b_process_band_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
BAND=process_band

# navigate to .../cluster_name
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# prevent implicit multithreading (otherwise tasks compete for resources)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# build the pipeline
pipetask build -p DRP-LoVoCCS.yaml#step3b \
    -C measure:${CLUSTER_DIR}/configs/coadd_3b_config_${BAND}.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step3b.yaml

# step3b measure
bps submit -b repo/repo \
    -i DECam/processing/quality_detectors_${BAND},DECam/processing/final_visit_summary_${BAND},DECam/processing/visit_summary_${BAND},DECam/processing/coadd_3a,refcats,skymaps \
    -o DECam/processing/coadd_3b_${BAND} \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step3b.yaml \
    -d "instrument='DECam' AND skymap='${CLN}_skymap' AND band='${BAND}'" \
    ${CLUSTER_DIR}/configs/bps_config_coadd3b.yaml

