#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH -J quality_check_cluster_name
#SBATCH -o slurm_outputs/quality_check_cluster_name-%j.out 
#SBATCH -e slurm_outputs/quality_check_cluster_name-%j.out 



#===========================
#AE: load the lsst science pipelines environment
LOAD_LSST="load_pipeline_path"
source ${LOAD_LSST}
setup lsst_distrib

#===========================
# Variables

script_folder="cluster_dir/python_scripts/quality_check"
output_folder="quality_check_output"
#cpunum="32"



#-----------------------

patches_tag="cluster_name_xminymin-xmaxymax"
#input_folder="photometric_correction_output"

color_term_residual_dered_csv="photometric_correction_output/${patches_tag}_match_star_refcat_dered_residual.csv"
stellar_locus_residual_dered_csv="photometric_correction_output/${patches_tag}_stellar_locus_dered_residual.csv"


#star_catalog0="photometric_correction_output/${patches_tag}_star_dered_dezp.csv"
#galaxy_catalog0="photometric_correction_output/${patches_tag}_gal_dered_dezp.csv"
all_catalog0="read_catalog_all_output/${patches_tag}_all.csv"
star_catalog0="photometric_correction_output/${patches_tag}_star.csv"
#star_catalog0="${output_folder}/${patches_tag}_star.csv"
galaxy_catalog0="photometric_correction_output/${patches_tag}_gal.csv"
before_clean_galaxy_catalog="photo_z_output/${patches_tag}_gal_dered_dezp_bpz_merge.csv"
cleaned_galaxy_catalog="${output_folder}/${patches_tag}_gal_clean.csv"


# We only consider mag0 shift here; no adjustment on magerr
all_catalog1="${output_folder}/${patches_tag}_all_dezp1.csv"
star_catalog1="${output_folder}/${patches_tag}_star_dezp1.csv"
galaxy_catalog1="${output_folder}/${patches_tag}_gal_dezp1.csv"
before_clean_star_catalog="${star_catalog1}"
cleaned_star_catalog="${output_folder}/${patches_tag}_star_clean.csv"


#selected_catalog=""
selected_all_catalog="${output_folder}/${patches_tag}_all_select.csv"
selected_star_catalog="${output_folder}/${patches_tag}_star_select.csv"
selected_galaxy_catalog="${output_folder}/${patches_tag}_gal_select.csv"
selected_cleaned_galaxy_catalog="${output_folder}/${patches_tag}_gal_clean_select.csv"
selected_cleaned_star_catalog="${output_folder}/${patches_tag}_star_clean_select.csv"


# Astrometry
#gaia_catalog="/gpfs/data/idellant/Clusters/ref_cat_database/gaia_dr2_cluster_name.csv"
gaia_catalog="/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new/cluster_name/gaia_dr3_cluster_name.csv"
#des_catalog="/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new/cluster_name/des_dr2_cluster_name.csv"
matched_star_gaia_catalog="${output_folder}/${patches_tag}_star_select_gaia.csv"
#matched_star_des_catalog="${output_folder}/${patches_tag}_star_select_des.csv"


# Correlation
correlation_selected_cleaned_galaxy_catalog="${output_folder}/${patches_tag}_gal_clean_select_corr.csv"
correlation_selected_cleaned_star_catalog="${output_folder}/${patches_tag}_star_clean_select_corr.csv"
#correlation_galaxy_galaxy_table="${output_folder}/${patches_tag}_gal_gal_corr.fits"
correlation_star_galaxy_table="${output_folder}/${patches_tag}_star_gal_corr.csv"
#correlation_star_star_table="${output_folder}/${patches_tag}_star_star_corr.fits"


# RS
spz_cat_name="photo_z_output/cluster_name_xminymin-xmaxymax_gal_match_specz_dered_dezp.csv"
pz_cat_name="photo_z_output/cluster_name_xminymin-xmaxymax_gal_dered_dezp_bpz_merge.csv"



#=======================
#-----------------------
# Script

[ ! -d ${output_folder} ] && mkdir ${output_folder} 


echo ""


####-----------------------
###echo "Star galaxy separation..."
###python ${script_folder}/separate_star_galaxy_catalog.py ${all_catalog0} ${star_catalog0} tmp
###
#-----------------------
# Note we only shift the zero point rather than the magerr.
# Because the goal is to check the efficiency of the detection.
# Note in photometric correction, we separated stars and galaxies. 
# Also for the galactic extinction, we don't need correct for it here. 
echo ""
echo "Zero point shift..."
python ${script_folder}/zero_point_correction2.py ${star_catalog0} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} ${star_catalog1}
python ${script_folder}/zero_point_correction2.py ${galaxy_catalog0} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} ${galaxy_catalog1}
python ${script_folder}/zero_point_correction2.py ${all_catalog0} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} ${all_catalog1}



#-----------------------
echo ""
echo "Selecting central sources for depth check..."

# Select sources near the field center
python ${script_folder}/select_central_region.py ${star_catalog1} ${selected_star_catalog} 
python ${script_folder}/select_central_region.py ${galaxy_catalog1} ${selected_galaxy_catalog} 

# This is for mag hist.
echo ""
echo "Selecting central sources for completeness check..."
python ${script_folder}/select_central_region.py ${all_catalog1} ${selected_all_catalog} 
####

#-----------------------
# This is for computing the correlation.
#echo "Cutting..."

echo ""
echo "Selecting sources for correlation check..."
echo ">>>>>>Cutting galaxy set..."
python ${script_folder}/cut2_galaxy.py ${before_clean_galaxy_catalog} ${cleaned_galaxy_catalog}

echo ">>>>>>Cutting star set..."
python ${script_folder}/cut2_star.py ${before_clean_star_catalog} ${cleaned_star_catalog}

###echo "Select central region..."
echo ">>>>>>>>>Selecting central galaxy set..."
python ${script_folder}/select_central_region.py ${cleaned_galaxy_catalog} ${selected_cleaned_galaxy_catalog} 

echo ">>>>>>>>>Selecting central star set..."
python ${script_folder}/select_central_region.py ${cleaned_star_catalog} ${selected_cleaned_star_catalog} 



#-----------------------
echo ""
echo "================================="
echo "Plotting mag vs mag err..."

python ${script_folder}/plot_mag_magerr.py ${selected_star_catalog} ${selected_galaxy_catalog}



#-----------------------

echo ""
echo "================================="
echo "Plotting hist of mag of fixed SNR & completeness..."
python ${script_folder}/plot_snr_mag_hist.py ${selected_star_catalog} ${selected_galaxy_catalog}

####python ${script_folder}/plot_mag_hist.py ${selected_star_catalog} ${selected_galaxy_catalog}
python ${script_folder}/plot_mag_hist.py ${selected_star_catalog} ${selected_galaxy_catalog} ${selected_all_catalog}



#-----------------------
echo ""
echo "================================="
echo "Checking the shape residual of PSF used stars..."
python ${script_folder}/check_star_shape.py ${all_catalog0} 



#-----------------------
echo ""
echo "================================="
echo "Computing correlation of shapes..."
# We use extra cuts on the catalogs before calculating the correlation
 
#echo "Further cutting the galaxy set..."
#python ${script_folder}/cut3_correlation.py ${selected_cleaned_galaxy_catalog} ${correlation_selected_cleaned_galaxy_catalog} 1
#echo "Further cutting the star set..."
#python ${script_folder}/cut3_correlation.py ${selected_cleaned_star_catalog} ${correlation_selected_cleaned_star_catalog} 0
# 
#
####echo "-----------------------------"
##echo "galaxy-galaxy correlation..."
##python ${script_folder}/correlation.py ${correlation_selected_cleaned_galaxy_catalog} ${correlation_selected_cleaned_galaxy_catalog} ${cpunum} ${correlation_galaxy_galaxy_table} 
##python ${script_folder}/plot_correlation.py ${correlation_galaxy_galaxy_table}
#
#
echo "Checking star-galaxy correlation..."
#python ${script_folder}/correlation.py ${correlation_selected_cleaned_star_catalog} ${correlation_selected_cleaned_galaxy_catalog} ${cpunum} ${correlation_star_galaxy_table} 
#python ${script_folder}/plot_correlation.py ${correlation_star_galaxy_table}
python ${script_folder}/correlation2.py ${selected_cleaned_star_catalog} ${selected_cleaned_galaxy_catalog} ${correlation_star_galaxy_table} 

#echo "star-star correlation..."
#python ${script_folder}/correlation.py ${correlation_selected_cleaned_star_catalog} ${correlation_selected_cleaned_star_catalog} ${cpunum} ${correlation_star_star_table} 
#python ${script_folder}/plot_correlation.py ${correlation_star_star_table}
#






#-----------------------
#-----------------------
echo ""
echo "================================="
echo "Comparing astrometry..."

python ${script_folder}/match_catalog.py "ra" "ra" "dec" "dec" "0.0002" "_cat" "_ref"  ${selected_star_catalog} ${gaia_catalog} ${matched_star_gaia_catalog} 
##python ${script_folder}/match_catalog.py "ra" "ra" "dec" "dec" "0.0002" "_cat" "_ref"  ${selected_star_catalog} ${des_catalog} ${matched_star_des_catalog} 

python ${script_folder}/compare_astrometry.py ${matched_star_gaia_catalog} cluster_name



#-----------------------
# Note if we want compare number counts, we need to constrain the sky region first.


#=======================
# Red sequence
# Consider r,i; g,r
echo ""
echo "================================="
# Find specz vs photometry cat, select member gal, list their photometry
# Find bpz vs photometry merge cat; then same
# Plot RS from final photometry cat
# Plot member gals with different marks

python ${script_folder}/red_sequence.py ${pz_cat_name} ${spz_cat_name}



#=======================
# Plot exposure table
python ${script_folder}/plot_exposure_table.py cluster_name

# Plot catalog
python ${script_folder}/plot_catalog_cmodel.py ${all_catalog0} 




#=======================
deactivate
