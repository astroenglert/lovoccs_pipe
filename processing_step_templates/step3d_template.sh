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
CLUSTER_DIR="/gpfs/data/idellant/englert_newPipelineDev/A85" # UPDATE LATER

# navigate to .../cluster_name
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# prevent implicit multithreading (otherwise tasks compete for resources)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# build the pipeline
pipetask build -p DRP-Merian_copy.yaml#step3d \
    -C selectGoodSeeingVisits:${CLUSTER_DIR}/configs/selectGoodSeeing_config.py \
    -C templateGen:${CLUSTER_DIR}/configs/templateGen_config.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step3d.yaml

# step3d create a template with high-quality seeing
# these exposures are passed from our first round of quality-checks...
# if this fails due to a lack of min-ccds, then we'll have to do some extra work to build the template from the calexps/visit summaries prior to check_visit
bps submit -b repo/repo \
    -i DECam/processing/coadd_3a,\
DECam/processing/final_visit_summary_g,DECam/processing/final_visit_summary_i,DECam/processing/final_visit_summary_r,DECam/processing/final_visit_summary_u,DECam/processing/final_visit_summary_z,\
skymaps \
    -o DECam/processing/coadd_3d \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step3d.yaml \
    -d "instrument='DECam' AND skymap='${CLN}_skymap'" \
    ${CLUSTER_DIR}/configs/bps_config.yaml


