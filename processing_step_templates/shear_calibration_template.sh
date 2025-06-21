#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=150GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -Jshear_calibration_cluster_name
#SBATCH -o slurm_outputs/shear_calibration_cluster_name-%j.out 
#SBATCH -e slurm_outputs/shear_calibration_cluster_name-%j.out 

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

# create an output directory for photometric_correction
mkdir shear_calibration_output


# this processing step just requires running calibrate_shears.py to produce the g1/g2
python -m python_scripts.shear_calibration.calibrate_shears "photo_z_output/${CLN}_dered_dezp_zphot_gals.csv" "shear_calibration_output/" "${CLN}_dered_dezp_zphot_scal_gals.csv" "naive"

