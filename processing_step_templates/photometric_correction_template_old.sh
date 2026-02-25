#!/bin/bash

#SBATCH --time=2:30:00
#SBATCH -n 1
#SBATCH --mem=15G
#SBATCH -J photometric_correction_cluster_name
#SBATCH -o slurm_outputs/photometric_correction_cluster_name-%j.out 
#SBATCH -e slurm_outputs/photometric_correction_cluster_name-%j.out 

#===========================
#AE: load the lsst science pipelines environment
LOAD_LSST="load_pipeline_path"
source ${LOAD_LSST}
setup lsst_distrib

#===========================
# Variables

#TODO update script directory
script_folder="cluster_dir/python_scripts/photometric_correction"
output_folder="photometric_correction_output"

#-----------------------

tag="all"
patches_tag="cluster_name_xminymin-xmaxymax"

ebv_image="${HOME}/data/Clusters/galactic_extinction_database/data_extinction_irsa/cluster_name.fits"
input_folder="read_catalog_all_output"
input_catalog_csv="${input_folder}/${patches_tag}_${tag}.csv"
dered_catalog_csv="${output_folder}/${patches_tag}_${tag}_dered.csv"
catalog_instrument="decam"

dered_catalog_star_csv="${output_folder}/${patches_tag}_star_dered.csv"
dered_catalog_gal_csv="${output_folder}/${patches_tag}_gal_dered.csv"
catalog_star_csv="${output_folder}/${patches_tag}_star.csv"
catalog_gal_csv="${output_folder}/${patches_tag}_gal.csv"

# For Skymapper we use dr2 because of its better photometry (but maybe fewer objects; in processing we might still use dr1)
#input_refcat_csv="/gpfs/data/idellant/Clusters/ref_cat_database/refcat_instrument_refcat_data_release_cluster_name.csv"

#AE: for gen3, we've relocated the refcats and python scripts, I've updated the filepaths accordingly
input_refcat_csv="/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new/cluster_name/refcat_instrument_refcat_data_release_cluster_name.csv"
dered_refcat_csv="${output_folder}/refcat_instrument_refcat_data_release_cluster_name_dered.csv"

matched_refcat_dered_catalog_csv="${output_folder}/${patches_tag}_match_star_refcat_dered.csv"


# To use color term, we need to match DECam observation and refcat first
# This file records residual mag and magerr
color_term_residual_dered_csv="${output_folder}/${patches_tag}_match_star_refcat_dered_residual.csv"

# To use DECam stellar locus, we can just look at the whole observation catalog (rather than the matched one)
# This file records residual mag and magerr
stellar_locus_residual_dered_csv="${output_folder}/${patches_tag}_stellar_locus_dered_residual.csv"


# Apply zero point correction
dered_catalog_gal_zero_point_corrected_csv="${output_folder}/${patches_tag}_gal_dered_dezp.csv"
dered_catalog_star_zero_point_corrected_csv="${output_folder}/${patches_tag}_star_dered_dezp.csv"
catalog_gal_zero_point_corrected_csv="${output_folder}/${patches_tag}_gal_dezp.csv"
catalog_star_zero_point_corrected_csv="${output_folder}/${patches_tag}_star_dezp.csv"


# Compare with DES
DES_catalog="/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new/cluster_name/des_dr1_cluster_name.csv"
matched_DES_dered_catalog_star_zero_point_corrected_csv="${output_folder}/${patches_tag}_star_match_DES_dered_dezp.csv"
matched_DES_catalog_star_zero_point_corrected_csv="${output_folder}/${patches_tag}_star_match_DES_dezp.csv"


# Compare with other catalogs



#=======================
#-----------------------
# Script

[ ! -d ${output_folder} ] && mkdir ${output_folder} 

echo "Patch: ${patches_tag}..."

echo ""
echo "Running apply_extinction_correction on input_catalog_csv: ${input_catalog_csv} ..."
# Note it applys to ALL objects and both PSF and CModel mags!
python ${script_folder}/apply_extinction_correction.py ${ebv_image} ${input_catalog_csv} ${dered_catalog_csv} ${catalog_instrument}



echo ""
echo "Running apply_extinction_correction on input_refcat_csv: ${input_refcat_csv}  ..."
python ${script_folder}/apply_extinction_correction.py ${ebv_image} ${input_refcat_csv} ${dered_refcat_csv} refcat_instrument



#=======================
# Here we divided the whole (dered) catalog into a star catalog and a galaxy catalog
python ${script_folder}/separate_star_galaxy_catalog.py ${dered_catalog_csv} ${dered_catalog_star_csv} ${dered_catalog_gal_csv}
# This is just for backup or quality check
python ${script_folder}/separate_star_galaxy_catalog.py ${input_catalog_csv} ${catalog_star_csv} ${catalog_gal_csv}



#-----------------------
echo ""
echo "================================="
echo "Matching catalog and refcat (stars)..."

# Note we can run a filter for stars first
echo "    Deredded catalog..."
if [ "refcat_instrument" == "sm" ];
    then
    #python ${script_folder}/match_catalog.py "ra" "raj2000" "dec" "dej2000" "0.0002" "_cat" "_ref"  ${dered_catalog_csv} ${dered_refcat_csv} ${matched_refcat_dered_catalog_csv} 
    python ${script_folder}/match_catalog.py "ra" "raj2000" "dec" "dej2000" "0.0002" "_cat" "_ref"  ${dered_catalog_star_csv} ${dered_refcat_csv} ${matched_refcat_dered_catalog_csv} 

elif [ "refcat_instrument" == "ps1" ];
    then
    python ${script_folder}/match_catalog.py "ra" "RAJ2000" "dec" "DEJ2000" "0.0002" "_cat" "_ref"  ${dered_catalog_star_csv} ${dered_refcat_csv} ${matched_refcat_dered_catalog_csv} 

else
    echo "Wrong refcat instrument!"
    exit 1
    
fi



#=======================



echo ""
echo "================================="
echo "Applying color-terms..."
echo ""
echo "*********Deredded catalog..."
# Note: to do this we need each object has both DECam and refcat mags
# The output offsets (and uncertainties) will be saved into another file (for BPZ)
# For convenience we don't input that to BPZ; we build a new catalog instead!
# Then in BPZ we only consider some CONSTANT and SMALL correction (which would be the same for most clusters?)

if [ "refcat_instrument" == "sm" ];
    then
    python ${script_folder}/apply_color_term_get_residual_sm.py ${matched_refcat_dered_catalog_csv} ${color_term_residual_dered_csv}
    
elif [ "refcat_instrument" == "ps1" ];
    then
    # Note for ps1 we can only "correct" grizY(y) bands
    python ${script_folder}/apply_color_term_get_residual_ps1.py ${matched_refcat_dered_catalog_csv} ${color_term_residual_dered_csv}



    echo ""
    echo "########################################"
    echo "Optional: For PS1, extra check on u-band, either use SDSS (deeper) or SkyMapper..."
    # We need extra check on u-band, either use sdss (deeper) or sm
    # Though we can skip this step, since we double check u-band by stellar locus (below; but could be less accurate)
    if [ -f "${HOME}/data/Clusters/calib_catalog_repo/catalogs_new/cluster_name/sdss_dr12_all_cluster_name.csv" ]; then
        input_refcat_csv="${HOME}/data/Clusters/calib_catalog_repo/catalogs_new/cluster_name/sdss_dr12_all_cluster_name.csv"
        dered_refcat_csv="${output_folder}/sdss_dr12_all_cluster_name_dered.csv"
        
        echo "--> apply_extinction_correction..."
        python ${script_folder}/apply_extinction_correction.py ${ebv_image} ${input_refcat_csv} ${dered_refcat_csv} sdss

        echo "--> match_catalog..."
        python ${script_folder}/match_catalog.py "ra" "RA_ICRS" "dec" "DE_ICRS" "0.0002" "_cat" "_ref"  ${dered_catalog_star_csv} ${dered_refcat_csv} ${matched_refcat_dered_catalog_csv}
        
        echo "--> apply_color_term_get_residual_sdss_u..."
        python ${script_folder}/apply_color_term_get_residual_sdss_u.py ${matched_refcat_dered_catalog_csv} ${color_term_residual_dered_csv}


    elif [ -f "${HOME}/data/Clusters/calib_catalog_repo/catalogs_new/cluster_name/sm_dr2_cluster_name.csv" ]; then 
        input_refcat_csv="${HOME}/data/Clusters/calib_catalog_repo/catalogs_new/cluster_name/sm_dr2_cluster_name.csv"
        dered_refcat_csv="${output_folder}/sm_dr2_cluster_name_dered.csv"
        
        echo "--> apply_extinction_correction..."
        python ${script_folder}/apply_extinction_correction.py ${ebv_image} ${input_refcat_csv} ${dered_refcat_csv} sm

        echo "--> match_catalog..."
        python ${script_folder}/match_catalog.py "ra" "raj2000" "dec" "dej2000" "0.0002" "_cat" "_ref"  ${dered_catalog_star_csv} ${dered_refcat_csv} ${matched_refcat_dered_catalog_csv}

        echo "--> apply_color_term_get_residual_sm_v..."
        python ${script_folder}/apply_color_term_get_residual_sm_v.py ${matched_refcat_dered_catalog_csv} ${color_term_residual_dered_csv}


    else
        echo "No good refcat for u-band! Skipping..."

    fi
        


else
    echo "Wrong refcat instrument!"
    exit 1
    
fi




echo ""
echo "================================="
echo "Getting offset from standard/model DECam stellar locus..."
# Correct catalogs with the biases from the above steps ("residual")
# Then calculate stellar locus (both data and std/theoretical colorterm curves)
# Especially u-band

#python ${script_folder}/compare_stellar_locus_decam_new.py ${dered_catalog_csv} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} 
python ${script_folder}/compare_stellar_locus_decam_new.py ${dered_catalog_star_csv} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} 



# Consider a final step: adjust mags, but not magerrs (because they are hard to add/subtract)
## Or leave this to photo-z
echo ""
echo "================================="
echo "Applying zero point correction..."
echo ""
echo "Running zero point correction on the original catalog..."
python ${script_folder}/zero_point_correction.py ${dered_catalog_gal_csv} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} ${dered_catalog_gal_zero_point_corrected_csv}
python ${script_folder}/zero_point_correction.py ${dered_catalog_star_csv} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} ${dered_catalog_star_zero_point_corrected_csv}
# This is just for quality check
python ${script_folder}/zero_point_correction.py ${catalog_gal_csv} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} ${catalog_gal_zero_point_corrected_csv}
python ${script_folder}/zero_point_correction.py ${catalog_star_csv} ${color_term_residual_dered_csv} ${stellar_locus_residual_dered_csv} ${catalog_star_zero_point_corrected_csv}



if [ -f "${DES_catalog}" ]; then

    echo ""
    echo "================================="
    echo "Comparing with DES..."
    echo ""
    # Compare with DES if possible
    python ${script_folder}/match_catalog.py "ra" "ra" "dec" "dec" 0.0002 "_dm" "_des" ${dered_catalog_star_zero_point_corrected_csv} ${DES_catalog} ${matched_DES_dered_catalog_star_zero_point_corrected_csv}
    python ${script_folder}/match_catalog.py "ra" "ra" "dec" "dec" 0.0002 "_dm" "_des" ${catalog_star_zero_point_corrected_csv} ${DES_catalog} ${matched_DES_catalog_star_zero_point_corrected_csv}
    
    # Plot difference
    echo ""
    echo "Checking DES: dered..."
    python ${script_folder}/plot_compare.py ${matched_DES_dered_catalog_star_zero_point_corrected_csv} 1

    echo ""
    echo "-----------------------------"
    echo "Checking DES: undered..."
    python ${script_folder}/plot_compare.py ${matched_DES_catalog_star_zero_point_corrected_csv} 0

fi


#-----------------------
# Added at Mon Jul 25 21:45:17 EDT 2022
# Here we consider using DECaLS (legacy_survey), PS1, SDSS, SM (we plot them all; their importance is sorted by this order) to plot the spatial distribution of magnitude difference.
# We use colorterms to convert the magnitudes, and remove a median magnitude offset. (similar to what we did in jointcal.sh to select the best order value) 
# Because the main goal is to check whether spatial pattern exists.

# Note these are multi-band catalogs becasue we need to convert their color to DECam magnitudes.
external_catalog_folder="/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new/cluster_name"
legacy_survey_catalog="${external_catalog_folder}/legacy_survey_dr9_cluster_name.csv"
ps1_catalog="${external_catalog_folder}/ps1_dr1_cluster_name.csv"
sdss_catalog="${external_catalog_folder}/sdss_dr12_all_cluster_name.csv"
sm_catalog="${external_catalog_folder}/sm_dr2_cluster_name.csv"

echo ""
echo ""
echo "================================="
echo "Comparing with other external catalogs..."
echo ""


#-----------------------
# DECaLS only has g, r, z bands
if [ -f ${legacy_survey_catalog} ]; then
    echo ""
    echo "---------------------------------------"
    echo "DECaLS catalog exists! Comparing..."
    
    external_catalog="${legacy_survey_catalog}"
    matched_catalog="${output_folder}/${patches_tag}_star_match_legacy_survey_dezp.csv"
    
    # Match catalog by RA/DEC
    # We use zp corrected DM catalog (for now) to how much difference compared with color-term converted catalog.
    # The difference should be large (if it's large then we need to find the reason).
    # Note a better option is to apply galactic extinction correction first (~0.01 mag difference).

    python ${script_folder}/match_catalog.py "ra" "ra" "dec" "dec" "0.00015" "_dm" "_legacy_survey" ${catalog_star_zero_point_corrected_csv} ${external_catalog} ${matched_catalog}

    for band in "g" "r" "z"; do

        # Plot difference
        python ${script_folder}/compare_mag_v3b.py ${matched_catalog} ${band} legacy_survey

    done

fi


#-----------------------
# PS1 has g,r,i,z,y bands
if [ -f ${ps1_catalog} ]; then
    echo ""
    echo "---------------------------------------"
    echo "PS1 catalog exists! Comparing..."
    
    external_catalog="${ps1_catalog}"
    matched_catalog="${output_folder}/${patches_tag}_star_match_ps1_dezp.csv"
    
    python ${script_folder}/match_catalog.py "ra" "RAJ2000" "dec" "DEJ2000" "0.00015" "_dm" "_ps1" ${catalog_star_zero_point_corrected_csv} ${external_catalog} ${matched_catalog}

    for band in "g" "r" "i" "z" "Y"; do

        # Plot difference
        python ${script_folder}/compare_mag_v3b.py ${matched_catalog} ${band} ps1

    done

fi


#-----------------------
# For SDSS we only consider u-band
if [ -f ${sdss_catalog} ]; then
    echo ""
    echo "---------------------------------------"
    echo "SDSS catalog exists! Comparing..."
    
    external_catalog="${sdss_catalog}"
    matched_catalog="${output_folder}/${patches_tag}_star_match_sdss_dezp.csv"
    
    python ${script_folder}/match_catalog.py "ra" "RA_ICRS" "dec" "DE_ICRS" "0.00015" "_dm" "_sdss" ${catalog_star_zero_point_corrected_csv} ${external_catalog} ${matched_catalog}

    for band in "u" "g" "r" "i" "z"; do

        # Plot difference
        python ${script_folder}/compare_mag_v3b.py ${matched_catalog} ${band} sdss

    done

fi


#-----------------------
# SM has v,g,r,i,z bands
if [ -f ${sm_catalog} ]; then
    echo ""
    echo "---------------------------------------"
    echo "SM catalog exists! Comparing..."
    
    external_catalog="${sm_catalog}"
    matched_catalog="${output_folder}/${patches_tag}_star_match_sm_dezp.csv"
    
    python ${script_folder}/match_catalog.py "ra" "raj2000" "dec" "dej2000" "0.00015" "_dm" "_sm" ${catalog_star_zero_point_corrected_csv} ${external_catalog} ${matched_catalog}

    for band in "u" "g" "r" "i" "z" "Y"; do

        # Plot difference
        python ${script_folder}/compare_mag_v3b.py ${matched_catalog} ${band} sm

    done

fi


#-----------------------
deactivate
