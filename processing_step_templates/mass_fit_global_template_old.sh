#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=10:00:00

# Default resources are 1 core with 2.8GB of memory.

# Use more memory (4GB):
#SBATCH --mem=96G

#SBATCH -n 20
#SBATCH -N 1

# Specify a job name:
#SBATCH -J mass_fit_global2_cluster_name

# Specify an output file
#SBATCH -o slurm_outputs/mass_fit_global_cluster_name-%j.out
#SBATCH -e slurm_outputs/mass_fit_global_cluster_name-%j.out

# Run a command

#=======================
# This script is based on the work of Jacqueline McCleary
# The idea is to fit all galaxies' shapes simultaneously with a mass profile (e.g. NFW)
# rather than (e.g. annular) binning first

# This script:
#   0. User need to provide a csv file for "peaks/centers"
#   1. Locate mass map peak
#   2. Fit mass with the peak coordinate

# The definitions of functions & constants here would be very similar to 
# (previous) annular-binning mass_fit script 
# (the structure of those functions & constants is learned from Nan Li)


#=======================
#AE: load the lsst science pipelines environment
LOAD_LSST="load_pipeline_path"
source ${LOAD_LSST}
setup lsst_distrib

#=======================
# Variables

cpun="96"

script_folder="cluster_dir/python_scripts/mass_fit_global"
output_folder="mass_fit_global_output"
#cluster_folder="/oscar/data/idellant/Clusters/cluster_name"
cluster_folder="."


#-----------------------
patches_tag="cluster_name_xminymin-xmaxymax"
input_folder="mass_map_output"

tag="gal"
#cleaned_catalog_csv="${input_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_cut.csv"
merged_catalog_g_csv="${input_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_cut_shear_calib_merge_cut.csv"
cluster_mass_map_peak_pixel_filename="${input_folder}/cluster_name_mass_map_peak_pixel.csv"




#-----------------------
# Just in case we would like adjust the peak file
tag1="_${RANDOM}"

#cluster_peak_radec_filename="./cluster_name_peak_radec.csv"
#cluster_peak_radec_pixel_filename="${output_folder}/cluster_name_peak_radec_pixel${tag1}.csv"
cluster_peak_pixel_filename="${output_folder}/cluster_name_peak_pixel${tag1}.csv"

# This is a list of masses from bootstrap
mass_filename="${output_folder}/${patches_tag}_${tag}_mass${tag1}.csv"
# Let the figure filename be mass_filename's ".csv" -> ".png"
mass_hist_figure="${output_folder}/${patches_tag}_${tag}_mass_hist${tag1}.png"


#-----------------------
# Script

[ ! -d ${output_folder} ] && mkdir ${output_folder} 



#-----------------------
#if [ ! -f ${cluster_peak_radec_filename} ]; then
#    echo "${cluster_peak_radec_filename} does not exist! Exiting..."
#    exit 1
#fi


#-----------------------
echo ""
echo "Checking peak radec file and building pixel file..."

# Note "cluster_peak_radec_pixel_filename" and "cluster_peak_pixel_filename" are different
# 1) One idea is to get peak(s) directly from mass map (maybe just the central & one peak)
# 2) Another idea is to get peak(s) from ra-dec coordinates and search for peaks on mass_map near these coordinates
# And then determine the x-y coordinates of the peaks
# 3) Or even simpler: directly convert ra-dec coordinates into x-y coordinates 
# 4) Or we can directly give x-y coordinates

#python ${script_folder}/radec2pixel.py ${cluster_folder} ${cluster_peak_radec_filename} ${cluster_peak_radec_pixel_filename}

# Leave this to future development
###python ${script_folder}/nearest_peak.py ${mass_map} ${cluster_peak_radec_pixel_filename} ${cluster_peak_pixel_filename}

#cp ${cluster_peak_radec_pixel_filename} ${cluster_peak_pixel_filename}
cp ${cluster_mass_map_peak_pixel_filename} ${cluster_peak_pixel_filename}




#-----------------------
echo ""
echo "Running mass_fit..."
python -u ${script_folder}/mass_fit_global_multi_center2_weight.py ${merged_catalog_g_csv} ${cluster_peak_pixel_filename} ${mass_filename} ${cpun} 
#python -u ${script_folder}/mass_fit_global_multi_center3.py ${merged_catalog_g_csv} ${cluster_peak_pixel_filename} ${mass_filename} ${cpun} 
#python -u ${script_folder}/mass_fit_global_multi_center3_landscape1D.py ${merged_catalog_g_csv} ${cluster_peak_pixel_filename} ${mass_filename} ${cpun} 



#-----------------------
echo ""
echo "hist.py..."
python ${script_folder}/hist_new_single.py ${mass_filename} ${mass_hist_figure}
# Added at Sat Jun 11 13:31:28 MST 2022
#python ${script_folder}/hist_new_1xN.py ${mass_filename} ${mass_hist_figure}



#-----------------------
deactivate
