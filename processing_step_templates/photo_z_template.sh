#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J photo_z_cluster_name
#SBATCH -o slurm_outputs/photo_z_cluster_name-%j.out 
#SBATCH -e slurm_outputs/photo_z_cluster_name-%j.out 

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

# load bash configs
source python_scripts/configs/processing_step_configs.sh

# create an output directory for photometric_correction
mkdir photo_z_output


# BPZ has two steps, first run bpz on a matched catalog w. spec-z's

# create a matched-catalog w. spec-z's
#TODO this needs to be overhauled to a standard-IO
#TODO here we cheat a little, specz ra/dec headers are the same as those used by default for decam (just ra/dec); should we create an "instrument" for ned-redshifts?
SPECZ_CSV="${SPECZ_DB}/${CLN}_ned_select.csv"

echo "Matching with the specz catalog"
python -m python_scripts.misc.match_catalogs "${SPECZ_CSV}" "decam" "_spec" "photometric_correction_output/${CLN}_dered_dezp_gals.csv" "decam" "" "photo_z_output/${CLN}_dered_dezp_gals_matched_specz.csv" "0.2"

# run bpz on the matched_catalog
echo "Running BPZ on the matched catalog"
python -m python_scripts.photo_z.bayesian_photo_z "photo_z_output/${CLN}_dered_dezp_gals_matched_specz.csv" "${CLN}_dered_dezp_gals_matched_specz_zphot.csv" "photo_z_output/" "20" "decam" "z_spec"


# Now run bpz on the entire catalog
echo "Running BPZ on the entire catalog"
python -m python_scripts.photo_z.bayesian_photo_z "photometric_correction_output/${CLN}_dered_dezp_gals.csv" "${CLN}_dered_dezp_zphot_gals.csv" "photo_z_output/" "20" "decam"



