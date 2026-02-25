#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J mass_map_cluster_name
#SBATCH -o slurm_outputs/mass_map_cluster_name-%j.out 
#SBATCH -e slurm_outputs/mass_map_cluster_name-%j.out 

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER

# prevent implicit multithreading (otherwise tasks compete for resources)
#export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# navigate to .../cluster_name
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# load filepaths from config
source python_scripts/configs/processing_step_configs.sh

# create an output directory for photometric_correction
mkdir mass_map_output


# this has two steps, run mass_map then draw the contours over the images!
#TODO eventually I'd like to find a way to remove the call to the coadd here
COADD="combine_patch_color_output/${CLN}_r33-88_deepCoadd.fits"
python -m python_scripts.mass_map.mass_map "shear_calibration_output/${CLN}_dered_dezp_zphot_scal_gals.csv" "100" "1.8" "mass_map_output/" "20" "decam" "1" "${CLN}" "schirmer"


# now draw contours
COADD_R="combine_patch_color_output/${CLN}_i33-88_deepCoadd.fits"
COADD_G="combine_patch_color_output/${CLN}_r33-88_deepCoadd.fits"
COADD_B="combine_patch_color_output/${CLN}_g33-88_deepCoadd.fits"
XRAY="${XRAY_DB}/${CLN}/Chandra/broad_flux_smoothed.fits"

# first just the mass_map over an inverse r-band image
python -m python_scripts.render_data.overlay_contours_map_inv "${COADD_R}" "mass_map_output/Map_E_peak_catalog_schirmer.csv" "mass_map_output/" "mass_map_output/${CLN}_r_inv_Map.png"

# then draw the mass_map over a color-image
python -m python_scripts.render_data.overlay_contours_map_rgb "${COADD_R}" "${COADD_G}" "${COADD_B}" "mass_map_output/Map_E_peak_catalog_schirmer.csv" "mass_map_output/" "mass_map_output/${CLN}_irg_Map.png"

# check for a chandra image and draw the x-ray contours if they exist
if [ -f "${XRAY}" ]; then
    python -m python_scripts.render_data.overlay_contours_map_inv "${COADD_R}" "mass_map_output/Map_E_peak_catalog_schirmer.csv" "mass_map_output/" "mass_map_output/${CLN}_r_inv_Map_xray.png" "${XRAY}"
fi


