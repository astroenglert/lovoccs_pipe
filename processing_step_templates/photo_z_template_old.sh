#!/bin/bash

#SBATCH --time=2:30:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=80G
#SBATCH -J photo_z_cluster_name
#SBATCH -o slurm_outputs/photo_z_cluster_name-%j.out 
#SBATCH -e slurm_outputs/photo_z_cluster_name-%j.out 



#===========================
# Get spec-z
# Run BPZ
# Compare photo-z with spec-z


#===========================
#AE: load the lsst science pipelines environment
LOAD_LSST="load_pipeline_path"
source ${LOAD_LSST}
setup lsst_distrib

#===========================
# Variables

#TODO update the script directory
cpu_num="20"
script_folder="cluster_dir/python_scripts/photo_z"
patches_tag="cluster_name_xminymin-xmaxymax"
input_folder="photometric_correction_output"
output_folder="photo_z_output"



#-----------------------

dered_catalog_gal_zero_point_corrected_csv="${input_folder}/${patches_tag}_gal_dered_dezp.csv"


specz_folder="/oscar/data/idellant/Clusters/spec_z_database"
specz_csv="${specz_folder}/cluster_name_ned_select.csv"

matched_specz_dered_catalog_zero_point_corrected_csv="${output_folder}/${patches_tag}_gal_match_specz_dered_dezp.csv"
matched_specz_dered_catalog_zero_point_corrected_bpz_csv="${output_folder}/${patches_tag}_gal_match_specz_dered_dezp_bpz.csv"


dered_catalog_gal_zero_point_corrected_bpz_csv="${output_folder}/${patches_tag}_gal_dered_dezp_bpz.csv"
dered_catalog_gal_bpz_merge_csv="${output_folder}/${patches_tag}_gal_dered_dezp_bpz_merge.csv"



#=======================
#-----------------------
# Script

[ ! -d ${output_folder} ] && mkdir ${output_folder} 


# Note in theory the contamination of stars in the catalog might affect the matching
# About selecting gals: maybe not useful since specz catalog only includes galaxies (including AGNs)
# and we can include star in the BPZ run and filter them out by extendedness later (but cost extra time)
echo ""
echo "Matching specz and catalog..."
python ${script_folder}/match_catalog.py "ra" "ra" "dec" "dec" "0.0002" "_cat" "_specz"  ${dered_catalog_gal_zero_point_corrected_csv} ${specz_csv} ${matched_specz_dered_catalog_zero_point_corrected_csv}



# Ingest the catalog into BPZ
# Note to run BPZ, we ONLY NEED to modify constant file (python script)
# So here we only modify constant file
echo ""
echo "Ingesting the matched catalog into BPZ..."
sed "s|data_filename|${matched_specz_dered_catalog_zero_point_corrected_csv}|g; s|output_filename|${matched_specz_dered_catalog_zero_point_corrected_bpz_csv}|g" ${script_folder}/constant_new2.py > constant_new.py



# Run BPZ
# Make link to the files and folders
# Run BPZ's shell script

echo ""
echo "Running BPZ..."

cp ${script_folder}/function_new.py function_new.py
cp ${script_folder}/prior_new.py prior_new.py

cp ${script_folder}/get_prob_new.py .
python -u get_prob_new.py ${cpu_num}

python ${script_folder}/bpz_specz.py ${matched_specz_dered_catalog_zero_point_corrected_bpz_csv} 



#=======================
# Run BPZ on the original catalog
echo ""
echo "Ingesting the original catalog into BPZ..."
sed "s|data_filename|${dered_catalog_gal_zero_point_corrected_csv}|g; s|output_filename|${dered_catalog_gal_zero_point_corrected_bpz_csv}|g" ${script_folder}/constant_new2.py > constant_new.py


#-----------------------
echo ""
echo "Running BPZ..."

cp ${script_folder}/get_prob_new.py .
python -u get_prob_new.py ${cpu_num}


#=======================
#-----------------------
# Merge the result with the original catalog (deredded)
echo ""
echo "Merging catalog and BPZ result..."

python ${script_folder}/horizontal_merge_table.py ${dered_catalog_gal_zero_point_corrected_csv} ${dered_catalog_gal_zero_point_corrected_bpz_csv} ${dered_catalog_gal_bpz_merge_csv} 


#-----------------------
# Check statistics of some terms
echo ""
echo "Plot stat of some terms in the merged (deredded original catalog + BPZ) catalog..."

python ${script_folder}/plot_stat.py ${dered_catalog_gal_bpz_merge_csv} 



#-----------------------
deactivate
