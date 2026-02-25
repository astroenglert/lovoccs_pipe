#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J red_sequence_cluster_name
#SBATCH -o slurm_outputs/red_sequence_cluster_name-%j.out 
#SBATCH -e slurm_outputs/red_sequence_cluster_name-%j.out 

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
mkdir red_sequence_output


# two steps here, run red_sequence.py then draw the contours
COADD="combine_patch_color_output/${CLN}_r33-88_deepCoadd.fits"
echo "mapping the rs..."
python -m python_scripts.red_sequence.red_sequence "photo_z_output/${CLN}_dered_dezp_zphot_gals.csv" "100" "1.8" "200" "red_sequence_output/" "decam" "${CLN}" "photo_z_output/${CLN}_dered_dezp_gals_matched_specz.csv" "z_spec"

# and finally render the density
echo "rendering contours..."
python -m python_scripts.render_data.overlay_contours_map_rs_inv "${COADD}" "red_sequence_output/rs_density_kde.fits" "mass_map_output/Map_E_peak_catalog_schirmer.csv" "mass_map_output/" "red_sequence_output/${CLN}_r_inv_Map_rs.png"

