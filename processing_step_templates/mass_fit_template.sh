#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J mass_fit_cluster_name
#SBATCH -o slurm_outputs/mass_fit_cluster_name-%j.out 
#SBATCH -e slurm_outputs/mass_fit_cluster_name-%j.out 

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

# create an output directory for photometric_correction
mkdir mass_fit_output

# this step just requires calling the corresponding python script
python -m python_scripts.mass_fit.mass_fit "shear_calibration_output/${CLN}_dered_dezp_zphot_scal_gals.csv" "mass_map_output/Map_E_peak_catalog_schirmer.csv" "${CLN}" "100" "mass_fit_output/" "20" "decam" "1" #(3,3) call to lower patch here is optional, but needs to be tweaked if different patches are used when computing the mass_map; default assumes 33-88... eventually we can drop this when mass_map is either running on-sky or outputting in a different wcs

