#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J meta_lensing_output_cluster_name
#SBATCH -o slurm_outputs/meta_lensing_output_cluster_name-%j.out
#SBATCH -e slurm_outputs/meta_lensing_output_cluster_name-%j.err

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

# create the output directories
mkdir meta_lensing_output
mkdir meta_lensing_output/calibrate
mkdir meta_lensing_output/mass_map
mkdir meta_lensing_output/mass_fit

# FIRST, run the calibrations!
python -m python_scripts.metadetect_shear.meta_calibrate "metadetect_processing/${CLN}_noshear_dered_dezp_zphot_gals.csv" "metadetect_processing/${CLN}_1p_dered_dezp_zphot_gals.csv" "metadetect_processing/${CLN}_1m_dered_dezp_zphot_gals.csv" "metadetect_processing/${CLN}_2p_dered_dezp_zphot_gals.csv" "metadetect_processing/${CLN}_2m_dered_dezp_zphot_gals.csv" "meta_lensing_output/calibrate/"


# NEXT, run mass_map!
COADD="combine_patch_color_output/${CLN}_r33-88_deepCoadd.fits"
python -m python_scripts.mass_map.mass_map "meta_lensing_output/calibrate/${CLN}_noshear_dered_dezp_zphot_gals_meta_cal.csv" "100" "1.8" "meta_lensing_output/mass_map/" "20" "decam" "1" "${CLN}" "schirmer"

# Also re-draw the corresponding contours
COADD_R="combine_patch_color_output/${CLN}_i33-88_deepCoadd.fits"
COADD_G="combine_patch_color_output/${CLN}_r33-88_deepCoadd.fits"
COADD_B="combine_patch_color_output/${CLN}_g33-88_deepCoadd.fits"
XRAY="${XRAY_DB}/${CLN}/Chandra/broad_flux_smoothed.fits"

# first just the mass_map over an inverse r-band image
python -m python_scripts.render_data.overlay_contours_map_inv "${COADD_R}" "meta_lensing_output/mass_map/Map_E_peak_catalog_schirmer.csv" "meta_lensing_output/mass_map/" "meta_lensing_output/mass_map/${CLN}_r_inv_Map.png"

# then draw the mass_map over a color-image
python -m python_scripts.render_data.overlay_contours_map_rgb "${COADD_R}" "${COADD_G}" "${COADD_B}" "meta_lensing_output/mass_map/Map_E_peak_catalog_schirmer.csv" "meta_lensing_output/mass_map/" "meta_lensing_output/mass_map/${CLN}_irg_Map.png"

# check for a chandra image and draw the x-ray contours if they exist
if [ -f "${XRAY}" ]; then
    python -m python_scripts.render_data.overlay_contours_map_inv "${COADD_R}" "meta_lensing_output/mass_map/Map_E_peak_catalog_schirmer.csv" "meta_lensing_output/mass_map/" "meta_lensing_output/mass_map/${CLN}_r_inv_Map_xray.png" "${XRAY}"
fi


# FINALLY, run mass_fit
python -m python_scripts.mass_fit.mass_fit "meta_lensing_output/calibrate/${CLN}_noshear_dered_dezp_zphot_gals_meta_cal.csv" "meta_lensing_output/mass_map/Map_E_peak_catalog_schirmer.csv" "${CLN}" "100" "meta_lensing_output/mass_fit/" "20" "decam" "0" #(3,3) call to lower patch here is optional, but needs to be tweaked if different patches are used when computing the mass_map; default assumes 33-88... eventually we can drop this when mass_map is either running on-sky or outputting in a different wcs


