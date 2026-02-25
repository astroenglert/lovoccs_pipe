#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=1:00:00

# Default resources are 1 core with 2.8GB of memory.

# Use more memory (4GB):
#SBATCH --mem=4G

# Specify a job name:
#SBATCH -J triaxiality_cluster_name

# Specify an output file
#SBATCH -o slurm_outputs/triaxiality_cluster_name-%j.out
#SBATCH -e slurm_outputs/triaxiality_cluster_name-%j.out

# Run a command


#=======================
# loading sextractor to get the bcg
module load intel-oneapi-compilers/2023.1.0-a7fw7qt
module load sextractor/2.25.0-wznqtm2

#AE: load the lsst science pipelines environment
LOAD_LSST="load_pipeline_path"
source ${LOAD_LSST}
setup lsst_distrib

#===========================
# Variables

#-----------------------
# FOLDERS

script_folder="cluster_dir/python_scripts/triaxiality"
output_folder="triaxiality_output"
#combine_patch_color_folder="combine_patch_color_output"
#mass_map_folder="mass_map_output"
#red_sequence_galaxy_distribution_folder="red_sequence_galaxy_distribution_output"


#-----------------------
# TAGS

patches_tag="cluster_name_xminymin-xmaxymax"


#-----------------------
# FILES

json_file="/oscar/data/idellant/Clusters/automatic_pipeline/red_sequence_galaxy_distribution/json_new/cluster_name.json"
#input_optical_image="${combine_patch_color_folder}/cluster_name_r55-66.fits"
#output_optical_image="${output_folder}/cluster_name_r55-66_cut.fits"
output_optical_image="${output_folder}/cluster_name_BCG_cut.fits"

BCG_record_filename="${output_folder}/BCG_record.csv"
BCG_angle_filename="${output_folder}/BCG_angle.csv"

#mass_map_peak_pixel_filename="${mass_map_folder}/cluster_name_mass_map_peak_pixel.csv"
mass_map_angle_filename="${output_folder}/mass_map_angle.csv"

RS_map_angle_filename="${output_folder}/RS_map_angle.csv"

r_cut_mpc="0.2"


#=======================
# We focus on the (2D) orientation of BCG, maps, galaxy distributions here.
# Note we don't measure elliptity here because it can vary with radial distance (also more difficult to measure).
echo ""
echo "Checking triaxiality..."

#-----------------------
# Script

[ ! -d ${output_folder} ] && mkdir ${output_folder} 


echo ""


#=======================
# BCG
# Since BCG can span different patch, and could be shredded in the detection, we run sextractor (either shell or python version) on coadd image (may cut near the center, use SIMBAD etc to get the BCG ra, dec).
# The image can be a small cut of the r-band coadd.
echo ""
echo "================================="
echo "Checking BCG..."


# 1. Make a cut of the r-band coadd near the NED cluster center.
#python ${script_folder}/crop_image_wcs2.py cluster_name ${input_optical_image} ${output_optical_image}
python ${script_folder}/crop_image_wcs2.py cluster_name ${json_file} ${output_optical_image} ${r_cut_mpc}
python ${script_folder}/fits_to_png.py ${output_optical_image}

# 2. Run sextractor on the cut.
echo "SNR,X,Y,THETA,AREA" > ${BCG_record_filename}




# 3. Match the sextractor catalog with NED catalog.
# This is to find BCG. Or maybe just the largest area one.

#for i in {5..100..5}; do
for i in {6..35..2}; do
#for i in {5..30..2}; do
#for i in {5..50..2}; do
#for i in {25..75..5}; do

    echo ""
    echo "-------------------------------------"
    echo "SNR=${i}..."

    sed "s/DET_TMP/${i}/g;s/PATH_TMP/${output_folder}\/test\.cat/g" ${script_folder}/tmp.config > ${output_folder}/use.config

    sex ${output_optical_image} -c ${output_folder}/use.config

#    ${script_folder}/sex2fiat ${output_folder}/test.cat > ${output_folder}/test.fiat

#    python ${script_folder}/fiat2csv2.py ${output_folder}/test.fiat ${output_folder}/test.csv

    python ${script_folder}/find_BCG.py ${output_folder}/test.cat ${i} >> ${BCG_record_filename}

done

python ${script_folder}/analyze_BCG_record.py ${BCG_record_filename} ${BCG_angle_filename} ${output_optical_image}


# 4. Turn the shape info to orientation.



#=======================
# Mass map
# We consider some typical SNR +- 10-20 percent difference.
# Compute orientation angle through second moments.
# The error bar is estimated by the change of the orientation angle when SNR changes.
# The typical value can be 50-80 percent of peak value.
# Need a radial cut e.g. < 0.5 Mpc for NED center (Xray should be the best; or maybe from mass map peak).
echo ""
echo "================================="
echo "Checking mass map..."

#-----------------------
# A possible filaments test: only look at things > 4-5 Mpc.
# Use 2nd moments or some better way to study that.

#Rs=$(python ${script_folder}/get_Rs.py ${mass_map_peak_pixel_filename})
#echo ""
#echo "Rs that maximizes peak SNR: ${Rs} pix......"
#mass_map_image="${mass_map_folder}/cluster_name_M_ap_SNR_b100_Rs${Rs}.fits"
#
#echo "Checking ${mass_map_image}..."

#python ${script_folder}/get_mass_map_angle.py ${mass_map_image} ${mass_map_angle_filename} cluster_name
python ${script_folder}/get_mass_map_angle.py ${json_file} ${mass_map_angle_filename} cluster_name 0

#-----------------------
# Another possible test is using quadrupole --- can we get orientation info from that?



#=======================
# RS galaxy distribution
# The first test is similar to mass map. We consider a typical value of some fraction of peak.
# Then we cut out the central region that have SNR +- some range. 
echo ""
echo "================================="
echo "Checking galaxy distribution..."

#density_image="${red_sequence_galaxy_distribution_folder}/cluster_name_density_b100_bw500.fits"

#echo "Checking ${density_image}..."

#python ${script_folder}/get_RS_map_angle.py ${density_image} ${RS_map_angle_filename} cluster_name
python ${script_folder}/get_RS_map_angle.py ${json_file} ${RS_map_angle_filename} cluster_name 1



#-----------------------
# Note another test can be using the data points of RS galaxies and directly compute second moments from points.
# Although the result should be quite similar to KDE's. Maybe noisier?



#=======================
deactivate
