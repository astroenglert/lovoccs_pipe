# source from the functions
source functions_Gen3.sh


# == RUNNING STEPS == #

# == STEP0 NOTES == # 
# After copying and pasting run_steps_Gen3 into a directory with the Cluster Name, 
# create_output creates a series of folders and python scripts. This only takes a few seconds to run.

# == STEP0 COMMAND == #
#create_output


# == STEP1 NOTES == # 
# This step queries and downloads frames from the NoirLab Science Archive.
# It usually takes a few hours to run, but if there is a lot of traffic or
# other bandwidth issues it can take as long as 12 hrs.

# == STEP1 COMMAND == #
#download_raw


# == STEP2 NOTES == #
# Unfortunately, LSSTPipe doesn't have any code built to catch errors due to
# missing bits in the fits header or other issues which may corrupt a fits file.
# As a result we have a step which manually checks each file (by opening each file with a python script
# and checking to see if the script crashes dramatically). This should only take a few minutes to run.

# == STEP2 COMMAND == #
#check_raws


# == STEP3 NOTES == # 
# Really... STEP2 is one script which I wrapped into a job-array to help things run a bit faster,
# but because of how it was written the corrupt fits have to be moved to a separate directory in another script.

# == STEP3 COMMAND == #
#move_corrupt_raws


# == STEP4 NOTES == #
# To begin LSSTPipe, we need to create a repository which will manage all of the files produced during processing.
# This is a quick process, so like the above steps this is run from the command line rather than with a batch script.

# == STEP4 COMMAND == #
#initialize_repo


# == STEP5 NOTES == #
# Next, we have to ingest our data into the repository. This step takes care of ingesting all of the catalogs,
# calibrations, and raw data into the repo. Depending on the amount of raw data, this can take an hour at most
# but usually takes less.

# == STEP5 COMMAND == #
#ingest_data


# == STEP6 NOTES == #
# This technically runs step1 of the Data Release Pipe (DRP), the big steps this includes are:
# - ISR (Instrumental Signal Removal: Crosstalk, Nonlinearity, Bias, Fringe, Flat, Brighter-Fatter)
# - Image Characterization (Initial Measurements, Cosmic-Ray Repair, Background Subtraction, Initial PSF Measurement)
# - Image Calibration (Initial Astrometry & Photometry)
# Remaining tasks consolidate sources into per-detector tables.
#
# For photometric refcats, the options are:
# ps1 (g,r,i,z,y), sdss (u-band only), sm_dr4 (g,i,r,u,z), see calib_catalog_repo for coverage
# This is the first BIG processing step, which takes a VERY long time to run, typically on the order of ~12hrs
# the argument formatting is CAT,BAND; e.g. to use ps1 for g,r; sm_dr4 for i; sm_dr4 for z; and sdss for u; the command is
# process_ccd sdss,u ps1,g ps1,r sm_dr4,i sm_dr4,z (the order of these is arbitrary; the command just runs a loop over the number of arguments)

# == STEP6 COMMAND == #
#process_ccd sm_dr4,u ps1,g ps1,r ps1,i ps1,z ps1,y


# == STEP7+8 NOTES == # 
# Before moving forward, we check the seeing in each visit/ccd and trim visits which are not "lensing-quality".
# We fit a Moffat profile to reference stars for measuring the FWHM and use the second moments to measure the distortion
# in the psf. Technically there is a function for doing this embedded in LSSTPipe, but we have yet to test it. After selecting
# the good detectors, a skymap object is created which encloses them.
# These usually take a few minutes to run each...

# == STEP7+8 COMMAND == #
#check_visit u g i r z
#select_visit


# == STEP9 NOTES == #
# This runs step2a of DRP. Taking the best CCD's following select_visit, we create final visit summaries
# which are required for jointcal. This only takes a few minutes to run...

# == STEP9 COMMAND == #
#visit_summary u g r i z


# == STEP10 NOTES == # 
# This runs step2b of DRP, joint-calibration. For precision astrophysics, errors in photometry/astrometry must be minimized.
# It happens that, on average, you can achive a more precise calibration by fitting the astrometry/photometry of each visit
# individually, rather than matching to a reference catalog after stacking. This can be done by fully modelling the
# atmosphere/optics (fgcm) or using an emiprical correction which models variations in brightnesses and positions
# between exposures (jointcal). Jointcal models variations in positions/brightnesses with a polynomial anchored in
# refcats and seeks to minimize a joint chi-squared. Ususally this takes about an hour-ish.
#
# Arguments should be IDENTICAL to the arguments used during check_visit, unless you are only trying to calibrate a specific band.

# == STEP10 COMMAND == #
#jointcal sm_dr4,u ps1,g ps1,r ps1,i ps1,z


# == STEP11 NOTES == # 
# This runs step2d of DRP. For reasons I can't quite remember, step2c is optional and isn't currently functional on DECam.
# It includes a final round of calibration which applies the corrections derived from jointcal and creates finalized
# visit summaries. This takes ~1hr/band with 20-cores.

# == STEP11 COMMAND == #
#final_visit_summary u g r i z


# == STEP12-15 NOTES == # 
# These steps run through step3a,step3b,step3c, and step3d of DRP.
# - step3a runs coaddition, detects sources, and runs the deblender (~3 hours w/ 140-cores)
# - step3b runs measurements on the coadds (~4 hrs/band w/ 20-cores)
# - step3c consolidates everything into tables and runs forced photometry (~5-7 hours w/ 140-cores)
# - step3d runs skycorrection, then makes a skycorrected coadded stack (~1 hour w/ 140-cores)
# Technically 3a,3b,3c can be run together if we disable the refcat matching... but since 3b is prone to errors
# it may not be very useful to merge them.

# == STEP12-15 COMMAND == #
#coadd_3a
#coadd_3b
#coadd_3c
#coadd_3d


# == STEP16 NOTES == #
# This step exports all of the data out of LSSTPipe, including:
# - Object catalogs containing shapes and magnitudes
# - Fits images containing the full fov
# - irg-images displaying the data
# Roughly ~1hr to run.

# == STEP16 COMMAND == #
#export_data


# == Processing with LSSTPipe is done and scripts only take a few minutes to run from here on == #


# == STEP17 NOTES == #
# From this point forward, the analysis is carried out on a catalog-level, the first step to this is
# correcting the photometry for color-terms and extinction (called de-redding). u-band is particularly funky here since
# there are very few existing refcats to calibrate the data with... so instead we use mock-observations
# created by integrating the spectra of reference-stars with the known transmission of the DECam filters
# and measure the u-band color terms directly from that calculation ("stellar-locus correction").
#
# Use "ps1" "dr2" or "sm" "dr4" depending on the primary refcat used during process_ccd.

# == STEP17 COMMAND == #
#photometric_correction "ps1" "dr2"


# == STEP18 NOTES == #
# Once we have a corrected catalog, it's time to identify galaxies in the field and estimate their redshift.
# For this we use a Bayesian PhotoZ algorithm.

# == STEP18 COMMAND == #
#photo_z


# == STEP19 NOTES == #
# By default, run shear calibration using HSC-Y1 calibration.
# In principle this is instrument dependent... but in practice HSC calib can be "good enough".
# We've implemented a modified version of metadetect (later processing steps) to provide a precise calibration.

# == STEP19 COMMAND == #
#shear_calibration


# == STEP20 NOTES == #
# Now that we have photo-z's, we can transform the HSM shape measurements into reduced shears and assemble a 
# mass_map. We use mass aperture statistics, which can be interpreted as:
# - An optimized matched-filter built to pick up NFW-like profiles
# - A measurement of the average convergence in an aperture, but weighted by a filter to maximize the signal

# == STEP20 COMMAND == #
#mass_map


# == STEP21 NOTES == #
# This step actually extracts the mass by fitting the shears of an NFW-profile to the observed shear-field and
# bootstrapping to produce error-bars.

# == STEP21 COMMAND == #
#mass_fit


# == STEP22 NOTES == #
# This stage looks for correlations between different objects and produces plots summarizing the quality of 
# observations/data which was used.

# == STEP22 COMMAND == #
#quality_check


# == STEP23 NOTES == #
# Generally, in a cluster the "red-sequence galaxies" are the oldest objects in the cluster and, as a result,
# effectively trace out the distribution of mass in the cluster. The red-sequence is selected from an HR-diagram
# of the cluster and smoothed contours representing their number density are produced.

# == STEP23 COMMAND == #
#red_sequence


# == Briefly going back to LSSTPipe for running metadetect! == #


# == STEP24-28 NOTES == #
# These steps run our implementation of metadetect.
# - meta_4a runs the shears on our coadds
# - meta_4b runs detect/deblend/measure on them, which we can use to build a robust calibration
# - meta_export collects all of the outputs from the 5 different versions of our coadds
# - meta_processing/meta_lensing finally process those coadds and run the lensing portion (including shear-calibration!)

# STEP24-28 COMMAND == #
#meta_4a
#meta_4b
#meta_export
#meta_processing
#meta_lensing


# == STEP29 NOTES == #
# Blast intermediate collections and other datasets, then zip the submit-directory.

# == STEP29 COMMAND == #
#gotta_blast


# And that's it!

# == WIP == #

#triaxiality "0" "0" "11" "11"
#subtract_stars

subtract_stars () {

	echo "Running STEP NA: subtract_stars"

	cp ${AUTO_PIPELINE_DIR}/config_templates/subtract_stars_process_config_template.py "${CLUSTER_DIR}/configs/subtract_stars_process_config.py"
	cp ${AUTO_PIPELINE_DIR}/config_templates/subtract_stars_measure_config_template.py "${CLUSTER_DIR}/configs/subtract_stars_measure_config.py"
	cp ${AUTO_PIPELINE_DIR}/config_templates/subtract_stars_subtract_config_template.py "${CLUSTER_DIR}/configs/subtract_stars_subtract_config.py"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/subtract_stars_template.sh > ${PROCESSING_STEP_DIR}/subtract_stars.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g" ${PROCESSING_STEP_DIR}/subtract_stars.sh
	
	# echo "Submitting to slurm..."
	# sbatch ${PROCESSING_STEP_DIR}/subtract_stars.sh
	# prompt_wait

}

# STEP NA: creating a template coadd with high-quality seeing for difference imaging!

coadd_3d () {

	echo "Running STEP 14: coadd_3d"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/step3d_template.sh > ${PROCESSING_STEP_DIR}/coadd_3d.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g" ${PROCESSING_STEP_DIR}/coadd_3d.sh
	
	# copying the config templates
	cp ${AUTO_PIPELINE_DIR}/config_templates/selectGoodSeeing_config_template.py "${CLUSTER_DIR}/configs/selectGoodSeeing_config.py"
	cp ${AUTO_PIPELINE_DIR}/config_templates/templateGen_config_template.py "${CLUSTER_DIR}/configs/templateGen_config.py"
	
	echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/coadd_3d.sh
	prompt_wait

}

# STEP NA: run the AP-Pipe to look for transients. 

ap_pipe () {

	echo "Running STEP 15: ap_pipe"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/ap_pipe_template.sh > ${PROCESSING_STEP_DIR}/ap_pipe.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g" ${PROCESSING_STEP_DIR}/ap_pipe.sh
	
	#INDEV might need to tweak some configs, so I'll leave these
	# copying the config templates
	#cp ${AUTO_PIPELINE_DIR}/config_templates/selectGoodSeeing_config_template.py "${CLUSTER_DIR}/configs/selectGoodSeeing_config.py"
	#cp ${AUTO_PIPELINE_DIR}/config_templates/templateGen_config_template.py "${CLUSTER_DIR}/configs/templateGen_config.py"
	
	echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/ap_pipe.sh
	prompt_wait

}


# STEP TODO: triaxiality
#TODO fix SExtractor configs so this works
triaxiality () {

    xmin=$1
    ymin=$2
    xmax=$3
    ymax=$4
    
    prompt "Running STEP 24: triaxiality"

    sed "s/cluster_name/${CLUSTER_NAME}/g; s/xmin/${xmin}/g; s/ymin/${ymin}/g; s/xmax/${xmax}/g; s/ymax/${ymax}/g" ${TEMPLATE_DIR}/triaxiality_template.sh > ${PROCESSING_STEP_DIR}/triaxiality.sh

    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g" ${PROCESSING_STEP_DIR}/triaxiality.sh

    echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/triaxiality.sh
    prompt_wait
    
    sleep 5s
 
}








