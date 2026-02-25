#!/bin/bash

#SBATCH --time=2:30:00
#SBATCH -n 1
#SBATCH --mem=15G
#SBATCH -J delta_sigma_profile_cluster_name
#SBATCH -o slurm_outputs/delta_sigma_profile_cluster_name-%j.out 
#SBATCH -e slurm_outputs/delta_sigma_profile_cluster_name-%j.out 



#===========================
#AE: load the lsst science pipelines environment
LOAD_LSST="load_pipeline_path"
source ${LOAD_LSST}
setup lsst_distrib

#===========================
# Variables

script_folder="cluster_dir/python_scripts/delta_sigma_profile"
output_folder="delta_sigma_profile_output"



#-----------------------

tag="gal"
patches_tag="cluster_name_xminymin-xmaxymax"

input_folder="mass_map_output"
#input_catalog_csv="${input_folder}/${patches_tag}_${tag}.csv"

merged_catalog_g_csv="${input_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_cut_shear_calib_merge_cut.csv"
cluster_mass_map_peak_pixel_filename="${input_folder}/cluster_name_mass_map_peak_pixel.csv"



#---------------------------
# Script

[ ! -d ${output_folder} ] && mkdir ${output_folder} 

echo "Making delta sigma profile..."
python ${script_folder}/delta_sigma_profile_xy.py ${merged_catalog_g_csv} ${cluster_mass_map_peak_pixel_filename} ${output_folder}
python ${script_folder}/delta_sigma_profile_xy_weight.py ${merged_catalog_g_csv} ${cluster_mass_map_peak_pixel_filename} ${output_folder}



#-----------------------
deactivate
