#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -J jointcal_process_band_cluster_name
#SBATCH -o slurm_outputs/jointcal_process_band_cluster_name-%j.out
#SBATCH -e slurm_outputs/jointcal_process_band_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PY_SCRIPTS="py_scripts"
BAND=process_band
PHOTOM=photom_ref

# navigate to .../cluster_name
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

# building the pipeline
pipetask build -p DRP-LoVoCCS.yaml#step2b \
    -C jointcal:${CLUSTER_DIR}/configs/jointcal_config_${BAND}.py \
    -s ${CLUSTER_DIR}/pipeline_yamls/DRP_step2b_${BAND}.yaml

# jointcal
# using pipetask didn't yield an increase in runtime AND caused read/write errors
# so going back to bps for this
bps submit -b repo/repo \
    -i DECam/calib/unbounded,DECam/processing/visit_summary_${BAND},refcats \
    -o DECam/processing/jointcal_${BAND} \
    -p ${CLUSTER_DIR}/pipeline_yamls/DRP_step2b_${BAND}.yaml \
    -d "instrument='DECam' AND band='${BAND}' AND skymap='${CLN}_skymap'" \
    ${CLUSTER_DIR}/configs/bps_config_jointcal.yaml

# jointcal residuals
if [ ${BAND} == "u" ] || [ ${BAND} == "g" ] || [ ${BAND} == "z" ]; then
    python python_scripts/render_data/jointcal_plot_residuals.py "${CLN}" "${PHOTOM}" 2 "${BAND}"
else
    python python_scripts/render_data/jointcal_plot_residuals.py "${CLN}" "${PHOTOM}" 1 "${BAND}"
fi

