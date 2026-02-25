#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J meta_processing_shear_type_cluster_name
#SBATCH -o slurm_outputs/meta_processing_shear_type_cluster_name-%j.out
#SBATCH -e slurm_outputs/meta_processing_shear_type_cluster_name-%j.err

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

# load in configs
source python_scripts/configs/processing_step_configs.sh

EXT_IMAGE="${EXT_DB}/${CLN}.fits"

# create an output
mkdir metadetect_processing

# loop over the five types of metadetect outputs
SHEARTYPE=shear_type

echo "Processing ${SHEARTYPE}"

INPUT_CAT="metadetect_export/${CLN}_metadetect_${SHEARTYPE}.csv"

# first run the extinction correction
echo "Running the extinction correction!"
python -m python_scripts.photometric_correction.extinction_correction ${EXT_IMAGE} ${INPUT_CAT} "metadetect_processing/${CLN}_${SHEARTYPE}_dered.csv" "decam"

# star-galaxy separation
echo "Separating stars and galaxies!"
python -m python_scripts.photometric_correction.separate_stars_galaxies "metadetect_processing/${CLN}_${SHEARTYPE}_dered.csv" "metadetect_processing/${CLN}_${SHEARTYPE}_dered_stars.csv" "metadetect_processing/${CLN}_${SHEARTYPE}_dered_gals.csv"

# applying zp-correction
echo "Applying ZP-correction to the galaxies!"
python -m python_scripts.photometric_correction.zero_point "metadetect_processing/${CLN}_${SHEARTYPE}_dered_gals.csv" "photometric_correction_output/${CLN}_matched_residuals.csv" "photometric_correction_output/${CLN}_matched_residuals_stellar_locus.csv" "metadetect_processing/${CLN}_${SHEARTYPE}_dered_dezp_gals.csv"

# Now run bpz on the catalog
echo "Running BPZ on the galaxies"
python -m python_scripts.photo_z.bayesian_photo_z "metadetect_processing/${CLN}_${SHEARTYPE}_dered_dezp_gals.csv" "${CLN}_${SHEARTYPE}_dered_dezp_zphot_gals.csv" "metadetect_processing/" "20" "decam"








