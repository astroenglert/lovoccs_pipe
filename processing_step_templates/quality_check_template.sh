#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J quality_check_cluster_name
#SBATCH -o slurm_outputs/quality_check_cluster_name-%j.out 
#SBATCH -e slurm_outputs/quality_check_cluster_name-%j.out 

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PY_SCRIPTS="py_scripts"
#TODO Pass this via run_steps rather than hard-coded
CALIB_CATALOG_REPO="/gpfs/data/idellant/Clusters/calib_catalog_repo" 

# navigate to .../cluster_name
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

# load configs
source python_scripts/configs/processing_step_configs.sh

# create an output directory for photometric_correction
mkdir quality_check_output


# quality check here is about the physical number of photons arriving at the detector, e.g. we care about quality of photometry after passing through the atmosphere/optics. Therefore, these steps don't use an extinction-corrected catalog.

# as a result there is a bit of work to be done; first separate stars and galaxies on the de-reddened catalog
#TODO if we change to a more sophisticated star-gal separation (e.g. check spectra rather than just extendedness), this will need heavy adjustments since it's not dereddened
INPUT_CAT="read_catalog_all_output/${CLN}_00-1111_all.csv"
echo "separating stars and galaxies"
python -m python_scripts.photometric_correction.separate_stars_galaxies "${INPUT_CAT}" "quality_check_output/${CLN}_stars.csv" "quality_check_output/${CLN}_gals.csv"


# next apply the zero-point correction to both catalogs
# really we apply the zp-correction here to be comparable w. refcats; but since we care about the signal post atm/optics, we skip adding the bias. I've created a modified version of zp-correction to handle this temporarily, really we should add a cln-option to the existing zp-correction script which disables the error-propagation

echo "running zp"
python -m python_scripts.photometric_correction.zero_point "quality_check_output/${CLN}_stars.csv" "photometric_correction_output/${CLN}_matched_residuals.csv" "photometric_correction_output/${CLN}_matched_residuals_stellar_locus.csv" "quality_check_output/${CLN}_dezp_stars.csv" "0"

python -m python_scripts.photometric_correction.zero_point "quality_check_output/${CLN}_gals.csv" "photometric_correction_output/${CLN}_matched_residuals.csv" "photometric_correction_output/${CLN}_matched_residuals_stellar_locus.csv" "quality_check_output/${CLN}_dezp_gals.csv" "0"


# next, match the star-catalog to Gaia for checking the photometry
REF_CAT="${CAT_DB}/${CLN}/gaia_dr3_${CLN}.csv"
echo "matching gaia"
python -m python_scripts.misc.match_catalogs "${REF_CAT}" "decam" "_ref" "quality_check_output/${CLN}_dezp_stars.csv" "decam" "_cat" "quality_check_output/${CLN}_dezp_stars_matched_gaia.csv" "0.2"


# now everything can be passed to quality_check, which takes care of the plotting in one big, somewhat messy, python script
echo "running quality_check.py"
python -m python_scripts.quality_check.quality_check "quality_check_output/${CLN}_dezp_stars.csv" "quality_check_output/${CLN}_dezp_gals.csv" "quality_check_output/${CLN}_dezp_stars_matched_gaia.csv" "_cat" "_ref" "query_result_${CLN}.csv" "quality_check_output/" "decam" "${CLN}"


