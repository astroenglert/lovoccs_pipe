#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J meta_4b_cluster_name
#SBATCH -o slurm_outputs/meta_4b_cluster_name-%j.out
#SBATCH -e slurm_outputs/meta_4b_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PY_SCRIPTS="py_scripts"
SHEARTYPE="shear_type"

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
# slightly lazy way of doing the configs... but its less work than writing 5x copies in the yaml
pipetask build -p DRP-LoVoCCS.yaml#step4b \
    -C measure:${CLUSTER_DIR}/configs/meta_4b_measure_config.py \
    -C forcedPhotCoadd:${CLUSTER_DIR}/configs/meta_4b_forced_config.py \
    -c detection:coaddName=deep_"${SHEARTYPE}"_ \
    -c detection:connections.inputCoaddName=deep_"${SHEARTYPE}"_ \
    -c detection:connections.outputCoaddName=deep_"${SHEARTYPE}"_ \
    -c detection:doScaleVariance=False \
    -c mergeDetections:coaddName=deep_"${SHEARTYPE}"_ \
    -c mergeDetections:connections.inputCoaddName=deep_"${SHEARTYPE}"_ \
    -c mergeDetections:connections.outputCoaddName=deep_"${SHEARTYPE}"_ \
    -c deblend:connections.inputCoaddName=deep_"${SHEARTYPE}"_ \
    -c deblend:connections.outputCoaddName=deep_"${SHEARTYPE}"_ \
    -c measure:coaddName=deep_"${SHEARTYPE}"_ \
    -c measure:connections.inputCoaddName=deep_"${SHEARTYPE}"_ \
    -c measure:connections.outputCoaddName=deep_"${SHEARTYPE}"_ \
    -c mergeMeasurements:coaddName=deep_"${SHEARTYPE}"_ \
    -c mergeMeasurements:connections.inputCoaddName=deep_"${SHEARTYPE}"_ \
    -c mergeMeasurements:connections.outputCoaddName=deep_"${SHEARTYPE}"_ \
    -c forcedPhotCoadd:coaddName=deep_"${SHEARTYPE}"_ \
    -c forcedPhotCoadd:connections.inputCoaddName=deep_"${SHEARTYPE}"_ \
    -c forcedPhotCoadd:connections.outputCoaddName=deep_"${SHEARTYPE}"_ \
    -c writeObjectTable:coaddName=deep_"${SHEARTYPE}"_ \
    -c writeObjectTable:connections.coaddName=deep_"${SHEARTYPE}"_ \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step4b_"${SHEARTYPE}".yaml

# only run metadetect on the central patches where data is deepest
PATCHES="(40..43, 51..56, 63..68, 75..80, 87..92, 100..103)"

# meta_4b
bps submit -b repo/repo \
    -i DECam/processing/meta_4a \
    -o DECam/processing/meta_4b_${SHEARTYPE} \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step4b_"${SHEARTYPE}".yaml \
    -d "instrument='DECam' AND skymap='${CLN}_skymap' AND (patch IN ${PATCHES})" \
    ${CLUSTER_DIR}/configs/bps_config.yaml

