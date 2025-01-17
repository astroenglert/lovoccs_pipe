#!/bin/bash 
# --- Start of slurm commands -----------

# Request an hour of runtime:
#SBATCH --time=11:59:59

# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=159G

# Specify a job name:
#SBATCH -J red_sequence_galaxy_distribution_cluster_name

# Specify an output file
# %j is a special variable that is replaced by the JobID when 
# job starts
#SBATCH -o red_sequence_galaxy_distribution_cluster_name-%j.out
#SBATCH -e red_sequence_galaxy_distribution_cluster_name-%j.out

#----- End of slurm commands ----

# Run a command


#===========================
#AE: load the lsst science pipelines environment
LOAD_LSST="load_pipeline_path"
source ${LOAD_LSST}
setup lsst_distrib

#===========================
# Variables

cpu_num="96"
script_folder="cluster_dir/python_scripts/red_sequence_galaxy_distribution"
output_folder="red_sequence_galaxy_distribution_output"
#DATA_path="/oscar/data/idellant/Clusters/cluster_name/DATA"
DATA_path="DATA"
#combine_patch_color_folder="/oscar/data/idellant/Clusters/cluster_name/combine_patch_color_output"
combine_patch_color_folder="combine_patch_color_output"
#SNR_threshold="3.5"



#-----------------------

patches_tag="cluster_name_xminymin-xmaxymax"
input_folder="photo_z_output"

tag="gal"
input_catalog_csv="${input_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge.csv"
gr_selected_catalog_csv="${output_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_RS_gr.csv"
ri_selected_catalog_csv="${output_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_RS_ri.csv"
gri_selected_catalog_csv="${output_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_RS_gri.csv"
final_selected_catalog_csv="${output_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_select.csv"
#shear_calibration_folder="/oscar/data/idellant/Clusters/automatic_pipeline/shear_calibration/hsc-y1-shear-calib"
#shear_calibration_output="${output_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_cut_shear_calib.csv"
#merged_catalog_csv="${output_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_cut_shear_calib_merge.csv"
#merged_catalog_g_csv="${output_folder}/${patches_tag}_${tag}_dered_dezp_bpz_merge_cut_shear_calib_merge_cut.csv"
cluster_red_sequence_galaxy_distribution_peak_pixel_filename="${output_folder}/cluster_name_red_sequence_galaxy_distribution_peak_pixel.csv"


#-----------------------
# "44-77" seems to be the best field size
patch_tag1="44-77_deepCoadd"
optical_image_R1="${combine_patch_color_folder}/cluster_name_i${patch_tag1}.fits" 
optical_image_G1="${combine_patch_color_folder}/cluster_name_r${patch_tag1}.fits" 
optical_image_B1="${combine_patch_color_folder}/cluster_name_g${patch_tag1}.fits" 

patch_tag2="33-88_deepCoadd"
optical_image_R2="${combine_patch_color_folder}/cluster_name_i${patch_tag2}.fits" 
optical_image_G2="${combine_patch_color_folder}/cluster_name_r${patch_tag2}.fits" 
optical_image_B2="${combine_patch_color_folder}/cluster_name_g${patch_tag2}.fits" 

Xray_image="/oscar/data/idellant/Clusters/Xray/cluster_name/Chandra/broad_flux_smoothed.fits"






#=======================
#-----------------------
# Script

[ ! -d ${output_folder} ] && mkdir ${output_folder} 



#-----------------------

echo ""
echo "Making selection on the catalog..."
python ${script_folder}/RS_select_gr.py ${input_catalog_csv} ${gr_selected_catalog_csv}
python ${script_folder}/RS_select_ri.py ${input_catalog_csv} ${ri_selected_catalog_csv}
python ${script_folder}/join_catalog.py ${gr_selected_catalog_csv} ${ri_selected_catalog_csv} idn ${gri_selected_catalog_csv}
python ${script_folder}/select_region.py ${gri_selected_catalog_csv} ${final_selected_catalog_csv}


##-----------------------
## Shear calibration
#echo ""
#echo "Running shear calibration..."
#
## Do shear calibration using some columns in the catalog.
#python ${shear_calibration_folder}/gen_hsc_calibrations.py ${selected_catalog_csv} ${shear_calibration_output} 
#
#
## Merge the shear calib result with the original catalog.
#python ${script_folder}/horizontal_merge_table.py ${selected_catalog_csv} ${shear_calibration_output} ${merged_catalog_csv}
#
#
## Make some figures for quality check, if needed.
#python ${script_folder}/cat_analysis.py ${merged_catalog_csv}
#
#
## 1) Get calibrated g; 2) Add a cut on the new g.
#python ${script_folder}/get_g.py ${merged_catalog_csv} ${merged_catalog_g_csv}




#-----------------------

echo ""
#echo "Building M_ap map and SNR map..."
echo "Building density map..."

python ${script_folder}/RS_KDE2.py ${final_selected_catalog_csv} 500 ${cpu_num} ${combine_patch_color_folder}
#python ${script_folder}/RS_KDE2.py ${final_selected_catalog_csv} 600 ${cpu_num} ${combine_patch_color_folder}
#python ${script_folder}/RS_KDE2.py ${final_selected_catalog_csv} 250 ${cpu_num} ${combine_patch_color_folder}

#for Rs in {3000..35000..1000}
#do
#    echo ""
#    echo "-----------------------------"
#    echo "--> Running Rs=${Rs}..."
#    python -W ignore ${script_folder}/schirmer_snr_weight.py ${merged_catalog_g_csv} ${Rs} ${cpu_num} ${DATA_path} ${SNR_threshold} 
##    python ${script_folder}/z_g_tomography.py ${merged_catalog_g_csv} ${output_folder}/cluster_name_M_ap_SNR_b100_Rs${Rs}_max.txt 
#done
#
#echo ""
#echo "Plotting Rs_vs_E_SNR_max..."
#python ${script_folder}/Rs_vs_E_SNR_max.py ${output_folder}/cluster_name_M_ap_SNR_b100 cluster_name ${cluster_red_sequence_galaxy_distribution_peak_pixel_filename}
#
#echo ""
#echo "Building cartoons..."
##convert -delay 100 -loop 0 ${output_folder}/cluster_name_M_ap_b100_Rs[1-9]*.png ${output_folder}/cluster_name_M_ap_b100_Rs.gif
#
#convert -delay 100 -loop 0 ${output_folder}/cluster_name_M_ap_SNR_b100_Rs[1-9]*000.png ${output_folder}/cluster_name_M_ap_SNR_b100_Rs.gif


#-----------------------
#echo ""
#echo "Building Significance map..."
#Rs="8000"
#python -W ignore ${script_folder}/schirmer_sig.py ${merged_catalog_g_csv} ${Rs} ${cpu_num}


#-----------------------
# Use Rs of max peak SNR to get the mass map fits image
# Then combine with optical images etc.
# Similar to combine fiatmap optical
#Rs=$(python ${script_folder}/get_Rs.py ${cluster_red_sequence_galaxy_distribution_peak_pixel_filename})
#echo ""
#echo "Rs that maximizes peak SNR: ${Rs} pix......"
#massmap_image_wcs="${output_folder}/cluster_name_M_ap_SNR_b100_Rs${Rs}.fits"
#massmap_image_wcs="mass_map_output/cluster_name_M_ap_SNR_b100_Rs12000.fits"
massmap_image_wcs="mass_map_output/cluster_name_M_ap_SNR_b100_Rs10000.fits"
density_image="${output_folder}/cluster_name_density_b100_bw500.fits"
#density_image="${output_folder}/cluster_name_density_b100_bw600.fits"
#density_image="${output_folder}/cluster_name_density_b50_bw500.fits"
#density_image="${output_folder}/cluster_name_density_b50_bw250.fits"

#-----------------------
echo "**********************"
echo "overplot_contour..."
#if [ -f ${Xray_image} ]; then
#    # Note if X-ray image exists, try this step:
#    # Note the cluster name is just for tag/title
#    python ${script_folder}/overplot_contour2.py ${optical_image_R} ${optical_image_G} ${optical_image_B} ${massmap_image_wcs} ${Xray_image} cluster_name ${patch_tag}
#    python ${script_folder}/overplot_contour2_inv.py ${optical_image_G} ${optical_image_G} ${optical_image_G} ${massmap_image_wcs} ${Xray_image} cluster_name ${patch_tag}
#
#else
#    python ${script_folder}/overplot_contour2.py ${optical_image_R} ${optical_image_G} ${optical_image_B} ${massmap_image_wcs} cluster_name ${patch_tag}
# Note image_G is actually r-band.
    python ${script_folder}/overplot_contour2_inv_test2.py ${optical_image_G1} ${optical_image_G1} ${optical_image_G1} ${massmap_image_wcs} ${density_image} cluster_name ${patch_tag1}
    python ${script_folder}/overplot_contour2_inv_test2.py ${optical_image_G2} ${optical_image_G2} ${optical_image_G2} ${massmap_image_wcs} ${density_image} cluster_name ${patch_tag2}

#fi



#echo ""
#
#echo "Montage..."
#python ${script_folder}/montage.py red_sequence_galaxy_distribution_output/cluster_name_M_ap_SNR_b100 ${SNR_threshold}

#-----------------------
deactivate
