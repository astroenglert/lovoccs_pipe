# Updating run_steps_5 to gen 3 LSP

# ==== #

# This script is built to work users through each task for processing clusters
# e.g. from raw -> calexp -> coadd -> catalogs -> lensing maps; each arrow represents a 'processing_step'
# Each step in the pipeline has a template file located at automatic_pipeline_v24/processing_step_templates [TEMP]
# run_steps copies each template file and uses text-replacement to pass the name of the cluster and other information into the template
# run_steps then submits the processing_step to slurm for exection
# The filled-in template are saved in data/$CLUSTER_NAME/processing_step

# == A quick note on the structure of outputs == #
# the previous version of run_steps dumped output files into the cluster directory, I've tried to ease that by creating a dedicated slurm_outputs folder. One small caviat of doing this is that, when a processing_step needs to be re-run, it must be sbatch'ed from the CLUSTER_NAME directory
# I've also added a feature which allocates resources, as necessary, with a little bit more intelligence. Check out https://ccv.brown.edu/rates/ to see how many cores and how much RAM you have available, and specify the default number of cores and RAM to be used for processing a cluster here!

CLUSTER_NAME='A85'
NODES=10 # might need to be tweaked for exploratory
CORES=180 # multiple of 10's only!
RAM=1300 # multiple of 10 only!
WALL_TIME=48:00:00 # default time is 2-days, may need to be adjusted for exploratory users

# == VARIABLES == #

#update when complete
AUTO_PIPELINE_DIR="/gpfs/data/idellant/Clusters/gen3_processing/testing_pipeline/A85_metadetect/lovoccs_pipe" #update when complete
TEMPLATE_DIR="${AUTO_PIPELINE_DIR}/processing_step_templates"
LOAD_PIPELINE_PATH="/gpfs/data/idellant/Clusters/gen3_processing/lsst_stack_v28_0_1/loadLSST.bash" # update to install in Clusters when complete
CLUSTER_DIR="/gpfs/data/idellant/Clusters/gen3_processing/testing_pipeline/A85_metadetect" #update when complete
PROCESSING_STEP_DIR="${CLUSTER_DIR}/processing_step"
CALIB_CATALOG_REPO="/gpfs/data/idellant/Clusters/calib_catalog_repo"


# == INITIALIZE LSP == #

# source ${LOAD_PIPELINE_PATH}
# setup lsst_distrib

# == FUNCTIONS ==  #

# helper function for customizing bps configs
bps_config_formatter () {

	## text-fill fields ##
	
	# photometry is passed to the formatter
	NUM_NODE="${1}"
	CORE_PER_NODE="${2}"
	MEM_PER_NODE="${3}"
	WALL_TIME="${4}"
	BPS_TAG="${5}"

	# running the actual text-replacement
	sed "s/NUM_NODE/${NUM_NODE}/g" ${AUTO_PIPELINE_DIR}/config_templates/bps_config_template.yaml > "${CLUSTER_DIR}/configs/bps_config${BPS_TAG}.yaml"
	sed -i "s/CORE_PER_NODE/${CORE_PER_NODE}/g" "${CLUSTER_DIR}/configs/bps_config${BPS_TAG}.yaml"
	sed -i "s/MEM_PER_NODE/${MEM_PER_NODE}/g" "${CLUSTER_DIR}/configs/bps_config${BPS_TAG}.yaml"
	sed -i "s/WALL_TIME/\"${WALL_TIME}\"/g" "${CLUSTER_DIR}/configs/bps_config${BPS_TAG}.yaml"

}

# helper function for printing text
prompt () {

	echo "-------------------------"
	echo $1
	echo "-------------------------"

}

# helper function to tell user a job has been submitted
prompt_wait () {

	prompt "Please wait for the sbatch job to finish! Use 'myq' and 'myjobinfo' to monitor progress."

}

# STEP 0: creates output directories for slurm-scripts and directory for processing_steps

create_output () {

echo "Running STEP 0: create_output"

# initialize .../CLUSTER_NAME
mkdir -p ${CLUSTER_DIR}/processing_step
mkdir -p ${CLUSTER_DIR}/slurm_outputs
mkdir -p ${CLUSTER_DIR}/configs
mkdir -p ${CLUSTER_DIR}/pipeline_yamls

cp -a ${AUTO_PIPELINE_DIR}/python_scripts/. ${CLUSTER_DIR}/python_scripts
cp ${AUTO_PIPELINE_DIR}/DRP-LoVoCCS.yaml ${CLUSTER_DIR}/DRP-LoVoCCS.yaml

echo "Setting up the default BPS-config"
echo "Using ${NODES} nodes, ${CORES} cores, and ${RAM} GB of RAM"
bps_config_formatter "${NODES}" "$((CORES/NODES))" "$((RAM/NODES))" "${WALL_TIME}"

}

# STEP 1: download raw data from NOAO

#TODO I've added a new script to manage the download using multiple streams in parallel
# BUT the connections are still throttled after a certain threshold is passed, we need to get on NoirLab to fix this!
download_raw () {

	echo "Running STEP 1: download_raw"

	echo "Downloading ${CLUSTER_NAME}..."
	echo "Raw images will be downloaded into .../${CLUSTER_NAME}/raw"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/noao_download_manager_template.sh > ${PROCESSING_STEP_DIR}/noao_download_manager.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/noao_download_manager.sh

	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/noao_download_manager.sh
	prompt_wait

}

# STEP 2: Often there are a few visits per cluster which have a bad fits header burried somewhere
# unfortunately LSSTPipe struggles a bit with these... so we have to check for them a-priori

check_raws () {

	echo "Running STEP 2: check_raws"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/check_raws_template.sh > ${PROCESSING_STEP_DIR}/check_raws.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/check_raws.sh

	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/check_raws.sh
	prompt_wait

}

# STEP 3: Now that we know which files are corrupt, we can move them to a separate folder
# in theory they may be repairable by fixing any bad bits or by at least recovering some of the CCD's
# but that's a project for someone in CS... not for me

move_corrupt_raws () {

	echo "Running STEP 3: move_corrupt_raws"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/move_corrupt_raws_template.sh > ${PROCESSING_STEP_DIR}/move_corrupt_raws.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/move_corrupt_raws.sh

	echo "Now running..."
	bash ${PROCESSING_STEP_DIR}/move_corrupt_raws.sh
	prompt_wait

}

# STEP 4: initialize the repository & create the Butler
# This, STEP 0, and move_corrupt_raws run directly where run_steps.sh is called
# usually this isn't a great thing to do since often it'll be run from a login node... 
# but neither step involve much processing, so it's safe to do so

initialize_repo () {

	echo "Running STEP 4: initialize_repo"
	cd "${CLUSTER_DIR}"
	
	echo "Setting up LSSTPipe"
	
	source ${LOAD_PIPELINE_PATH}
	setup lsst_distrib
	
	echo "Creating the repo and setting up DECam"
	
	# tell LSSTPipe to create a repository here and register it for DECam data
	butler create repo/repo
	butler register-instrument repo/repo lsst.obs.decam.DarkEnergyCamera
	
	# write instrument-specific calibration frames to the repository (only a few MB)
	butler write-curated-calibrations --label curated repo/repo lsst.obs.decam.DarkEnergyCamera

}

# STEP 5: ingest raws, calibrations, and catalogs to the repository

ingest_data () {

	echo "Running STEP 5: ingest_data"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/ingest_data_template.sh > ${PROCESSING_STEP_DIR}/ingest_data.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/ingest_data.sh

	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/ingest_data.sh
	prompt_wait

}

# STEP 6: perform instrumental signal removal (flats/bias/xtalk) and perform a rough astrom/photom calibration

# config-formatting helper function

process_ccd_refcat_formatter () {
	
	## text-fill fields ##
	
	# gaia is hard-coded as the default for astrometry
	ASTROM_REF="gaia"
	ASTROM_BAND="phot_g_mean_mag"
	
	# photometry is passed to the formatter
	PHOTOM_REF="${1}"
	PHOTOM_MAG="${2}"
	
	# last argument is the physical band being processed
	PROCESS_BAND="${3}"

	# running the actual text-replacement
	sed "s/astrom_ref/${ASTROM_REF}/g" ${AUTO_PIPELINE_DIR}/config_templates/process_ccd_refcat_config_template.py > "${CLUSTER_DIR}/configs/process_ccd_refcat_config_${PROCESS_BAND}.py"
	sed -i "s/astrom_band/${ASTROM_BAND}/g" "${CLUSTER_DIR}/configs/process_ccd_refcat_config_${PROCESS_BAND}.py"
	sed -i "s/photom_ref/${PHOTOM_REF}/g" "${CLUSTER_DIR}/configs/process_ccd_refcat_config_${PROCESS_BAND}.py"
	sed -i "s/photom_mag/${PHOTOM_MAG}/g" "${CLUSTER_DIR}/configs/process_ccd_refcat_config_${PROCESS_BAND}.py"
	sed -i "s/process_band/${PROCESS_BAND}/g" "${CLUSTER_DIR}/configs/process_ccd_refcat_config_${PROCESS_BAND}.py"
	
}

process_ccd_psf_formatter () {
	
	# the physical band being processed
	PROCESS_BAND="${1}"
	U_FLUX_LIM="False"
	U_STD_WIDTH=1
	U_PSF_ITER=1
	U_MIN_FLUX=100
	DEFAULT_FLUX=4000
	PSF_MODEL='psfex'
	
	# pass different min-flux for psf-candidates depending on the band
	# u-band is generally much fainter than rgizY
	
	if [ ${PROCESS_BAND} == 'u' ]; then
		sed "s/min_flux/${U_MIN_FLUX}/g" ${AUTO_PIPELINE_DIR}/config_templates/process_ccd_psf_config_template.py > "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
		sed -i "s/std_width/${U_STD_WIDTH}/g" "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
		sed -i "s/psf_iter/${U_PSF_ITER}/g" "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
		sed -i "s/flux_lim/${U_FLUX_LIM}/g" "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
	else
		sed "s/min_flux/${DEFAULT_FLUX}/g" ${AUTO_PIPELINE_DIR}/config_templates/process_ccd_psf_config_template.py > "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
		sed -i "s/std_width/0.15/g" "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
		sed -i "s/psf_iter/2/g" "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
		sed -i "s/flux_lim/True/g" "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
	fi

	# update the psf-model
	sed -i "s/psf_determiner/${PSF_MODEL}/g" "${CLUSTER_DIR}/configs/process_ccd_psf_config_${PROCESS_BAND}.py"
	
}

process_ccd () {

	echo "Running STEP 6: process_ccd"
	
	# because different clusters pull from various catalogs, separate configs are created for each band
	# the for-loops here are a little lazy... but there isn't any tangible benefit to handling it more neatly
	# CATALOG is the photometry catalog, BAND is the band to process using this catalog
	for INPUT in "$@"; do
		CATALOG=${INPUT/%,*/}
		BAND=${INPUT/#*,/}
		
		# ps1 formatting
		if [ "${CATALOG}" == "ps1" ]; then
			PS1_BANDS=('g' 'r' 'i' 'z' 'y')
			PS1_MAPS=('gmag' 'rmag' 'imag' 'zmag' 'ymag')
			
			for i in ${!PS1_BANDS[@]}; do
				if [ "${PS1_BANDS[$i]}" == "$BAND" ]; then 
					process_ccd_refcat_formatter "${CATALOG}" "${PS1_MAPS[$i]}" "${BAND}"
					process_ccd_psf_formatter "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: PS1 does not contain ${BAND}"; return; fi
			done
		fi
		
		# sky-mapper dr1 formatting
		if [ "${CATALOG}" == "sm_dr1" ]; then
			SM1_BANDS=('u' 'g' 'r' 'i' 'z')
			SM1_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM1_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM1_BANDS[@]}; do
				if [ "${SM1_BANDS[$i]}" == "$BAND" ]; then 
					process_ccd_refcat_formatter "${CATALOG}_${SM1_CAT_BANDS[$i]}" "${SM1_MAPS[$i]}" "${BAND}"
					process_ccd_psf_formatter "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR1 does not contain ${BAND}"; return; fi
			done
		fi
		
		# sky-mapper dr2 formatting
		if [ "${CATALOG}" == "sm_dr2" ]; then
			SM2_BANDS=('u' 'g' 'r' 'i' 'z')
			SM2_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM2_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM2_BANDS[@]}; do
				if [ "${SM2_BANDS[$i]}" == "$BAND" ]; then 
					process_ccd_refcat_formatter "${CATALOG}_${SM2_CAT_BANDS[$i]}" "${SM2_MAPS[$i]}" "${BAND}"
					process_ccd_psf_formatter "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR2 does not contain ${BAND}"; return; fi
			done
		fi
		
		# sdss formatting
		if [ "${CATALOG}" == "sdss" ]; then
			SDSS_BANDS=('u')
			SDSS_MAPS=('upmag')
			
			for i in ${!SDSS_BANDS[@]}; do
				if [ "${SDSS_BANDS[$i]}" == "$BAND" ]; then 
					process_ccd_refcat_formatter "${CATALOG}" "${SDSS_MAPS[$i]}" "${BAND}"
					process_ccd_psf_formatter "${BAND}"
					break
				fi
				if [ $i == 0 ]; then echo "ERROR: SDSS does not contain ${BAND}"; return; fi
			done
		fi
		
		if [ "${CATALOG}" == "des_dr2" ]; then
			DES_BANDS=('g' 'r' 'i' 'z' 'y')
			DES_MAPS=('wavg_mag_psf_g' 'wavg_mag_psf_r' 'wavg_mag_psf_i' 'wavg_mag_psf_z' 'wavg_mag_psf_y')
			
			for i in ${!DES_BANDS[@]}; do
				if [ "${DES_BANDS[$i]}" == "$BAND" ]; then 
					process_ccd_refcat_formatter "${CATALOG}" "${DES_MAPS[$i]}" "${BAND}"
					process_ccd_psf_formatter "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: DES does not contain ${BAND}"; return; fi
			done
		fi

		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/process_ccd_template.sh > ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh

		# echo "Submitting to slurm..."
		# sbatch ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh
		# sleep 15m

		# Soren's idea, submit successive jobs so that they wait in the queue for the previous to finish.
		# This let's us better utilize resources, at the cost of some overhead ran in serial at the beginning/end of each job
		if [ -z "$JOBID" ]; then
			JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh)
			echo "Submitted process_ccd_${BAND} with ${JOBID}"
		else
			JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh)
			echo "Submitted process_ccd_${BAND} with ${JOBID}"
		fi
	done

}

# STEP 7: check each visit to measure the seeing (PSF-FWHM) and average ellipticity across the sensor

# TODO: Test the seeing-check in the LSP against our algorithm (Moffat-fits)
# Currently in the LSP, there is a step built into the coaddition stage which will pick-out good seeing visits
# How well does this step work and how does it compare to our method of selecting detectors/visits with good seeing?

check_visit () {
	
	echo "Running STEP 7: check_visit"

	FWHM=6.4
	ELLIP=0.33
	FWHM_r=4.4
	ELLIP_r=0.13

	for BAND in "$@"; do
		
		case ${BAND} in
		"u" | "g" | "i" | "z" | "y")
		
		# pass the cluster name into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/check_visit_template.sh > ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		# pass the current lsst_pipeline into the template
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		# pass the band being processed into the template
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		sed -i "s/fwhm_cut/${FWHM}/g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh		
		sed -i "s/ellip_cut/${ELLIP}/g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
				
		echo "Submitting check_visit_${BAND} slurm..."
		sbatch ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		
		;;
		"r")
		
		# pass the cluster name into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/check_visit_template.sh > ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		# pass the current lsst_pipeline into the template
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		# pass the band being processed into the template
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		sed -i "s/fwhm_cut/${FWHM_r}/g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh		
		sed -i "s/ellip_cut/${ELLIP_r}/g" ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
				
		echo "Submitting check_visit_${BAND} slurm..."
		sbatch ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		
		;;
		*)

		echo "${BAND} is not allowed... please choose between ugirzy!"
		;;
		esac
	done

	prompt_wait
}

# STEP 8: apply star-count/fwhm/ellip cuts and create a skymap enclosing our exposures
# this is the second half of "check_visit", which is separate since it cannot be parallelized

select_visit () {

	echo "Running STEP 8: select_visit"

	# fwhm/ellip-cuts are band dependent
	FWHM=6.4
	ELLIP=0.33
	FWHM_r=4.4
	ELLIP_r=0.13

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/select_visit_template.sh > ${PROCESSING_STEP_DIR}/select_visit.sh

	# pass the current lsst_pipeline into the template
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/select_visit.sh

	# passing the fwhm/ellip cuts to select_visit
	sed -i "s/fwhm_cut_r/${FWHM_r}/g" ${PROCESSING_STEP_DIR}/select_visit.sh
	sed -i "s/ellip_cut_r/${ELLIP_r}/g" ${PROCESSING_STEP_DIR}/select_visit.sh
	sed -i "s/fwhm_cut/${FWHM}/g" ${PROCESSING_STEP_DIR}/select_visit.sh
	sed -i "s/ellip_cut/${ELLIP}/g" ${PROCESSING_STEP_DIR}/select_visit.sh

	# copying the skymap configs
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${AUTO_PIPELINE_DIR}/config_templates/skymap_config_template.py > ${CLUSTER_DIR}/configs/skymap_config.py

	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/select_visit.sh
	prompt_wait

}

# STEP 9: consolidate tables from process_ccd and ready data for jointcal

visit_summary () {

	echo "Running STEP 9: visit_summary"
	
	for BAND in "$@"; do
		
		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/visit_summary_template.sh > ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh

		# echo "Submitting to slurm..."
		# sbatch ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		# sleep 1m
		
		# Soren's idea, submit successive jobs so that they wait in the queue for the previous to finish.
		# This let's us better utilize resources, at the cost of some overhead ran in serial at the beginning/end of each job
		if [ -z "$JOBID" ]; then
			JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh)
			echo "Submitted visit_summary_${BAND} with ${JOBID}"
		else
			JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh)
			echo "Submitted visit_summary_${BAND} with ${JOBID}"
		fi


	done

}

# STEP #10 joint-calibration, compares visits of the same band to refcat and carries out a joint photom/astrom-cal

# config-formatting helper function

jointcal_formatter () {
	
	## text-fill fields ##
	
	# gaia is hard-coded as the default for astrometry
	ASTROM_REF="gaia"
	ASTROM_BAND="phot_g_mean_mag"
	
	# photometry is passed to the formatter
	PHOTOM_REF="${1}"
	PHOTOM_MAG="${2}"
	
	# last argument is the physical band being processed
	PROCESS_BAND="${3}"
	
	# setting the order of jointcal based on the band, easy to tweak as necessary

	if [ "${PROCESS_BAND}" == "u" ] || [ "${PROCESS_BAND}" == "g" ] || [ "${PROCESS_BAND}" == "z" ]; then
		JOINTCAL_ORDER=2 # ugz have best results at order 2 (at least for A85...)
	else 
		JOINTCAL_ORDER=1 # ri have best results at order 1... need to test this!
	fi
	
	# running the actual text-replacement
	sed "s/astrom_ref/${ASTROM_REF}/g" ${AUTO_PIPELINE_DIR}/config_templates/jointcal_config_template.py > "${CLUSTER_DIR}/configs/jointcal_config_${PROCESS_BAND}.py"
	sed -i "s/astrom_band/${ASTROM_BAND}/g" "${CLUSTER_DIR}/configs/jointcal_config_${PROCESS_BAND}.py"
	sed -i "s/photom_ref/${PHOTOM_REF}/g" "${CLUSTER_DIR}/configs/jointcal_config_${PROCESS_BAND}.py"
	sed -i "s/photom_mag/${PHOTOM_MAG}/g" "${CLUSTER_DIR}/configs/jointcal_config_${PROCESS_BAND}.py"
	sed -i "s/process_band/${PROCESS_BAND}/g" "${CLUSTER_DIR}/configs/jointcal_config_${PROCESS_BAND}.py"
	sed -i "s/jointcal_order/${JOINTCAL_ORDER}/g" "${CLUSTER_DIR}/configs/jointcal_config_${PROCESS_BAND}.py"
}

jointcal () {

	echo "Running STEP 10: jointcal"
	
	# let's get the new configs for bps in CLN
	cp ${AUTO_PIPELINE_DIR}/config_templates/bps_config_jointcal.yaml "${CLUSTER_DIR}/configs/bps_config_jointcal.yaml"
	
	# separate directory for debug-data (residuals)
	for i in {0..7}; do
            mkdir -p "jointcal_residuals/order_${i}"
	done
	
	# we use a custom bps_config to let tasks run in parallel
	# unfortunately jointcal is a serial process, so this can't be spread up by much
	echo "Setting up the jointcal BPS-config"
	echo "Using 1 node, 10 cores, and 90 GB of RAM"
	bps_config_formatter "1" "10" "90" "48:00:00" "_jointcal"
	
	# because different clusters pull from various catalogs, separate configs are created for each band
	# the for-loops here are a little lazy... but there isn't any tangible benefit to handling it more neatly
	# CATALOG is the photometry catalog, BAND is the band to process using this catalog
	for INPUT in "$@"; do
		CATALOG=${INPUT/%,*/}
		BAND=${INPUT/#*,/}
		
		# ps1 formatting
		if [ "${CATALOG}" == "ps1" ]; then
			PS1_BANDS=('g' 'r' 'i' 'z' 'y')
			PS1_MAPS=('gmag' 'rmag' 'imag' 'zmag' 'ymag')
			
			for i in ${!PS1_BANDS[@]}; do
				if [ "${PS1_BANDS[$i]}" == "$BAND" ]; then 
					jointcal_formatter "${CATALOG}" "${PS1_MAPS[$i]}" "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: PS1 does not contain ${BAND}"; return; fi
			done
		fi
		
		# sky-mapper dr1 formatting
		
		
		
		if [ "${CATALOG}" == "sm_dr1" ]; then
			SM1_BANDS=('u' 'g' 'r' 'i' 'z')
			SM1_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM1_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM1_BANDS[@]}; do
				if [ "${SM1_BANDS[$i]}" == "$BAND" ]; then 
					jointcal_formatter "${CATALOG}_${SM1_CAT_BANDS[$i]}" "${SM1_MAPS[$i]}" "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR1 does not contain ${BAND}"; return; fi
			done
		fi
		
		# sky-mapper dr2 formatting
		if [ "${CATALOG}" == "sm_dr2" ]; then
			SM2_BANDS=('u' 'g' 'r' 'i' 'z')
			SM2_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM2_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM2_BANDS[@]}; do
				if [ "${SM2_BANDS[$i]}" == "$BAND" ]; then 
					jointcal_formatter "${CATALOG}_${SM2_CAT_BANDS[$i]}" "${SM2_MAPS[$i]}" "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR2 does not contain ${BAND}"; return; fi
			done
		fi
		
		# sdss formatting
		if [ "${CATALOG}" == "sdss" ]; then
			SDSS_BANDS=('u')
			SDSS_MAPS=('upmag')
			
			for i in ${!SDSS_BANDS[@]}; do
				if [ "${SDSS_BANDS[$i]}" == "$BAND" ]; then 
					jointcal_formatter "${CATALOG}" "${SDSS_MAPS[$i]}" "${BAND}"
					break
				fi
				if [ $i == 0 ]; then echo "ERROR: SDSS does not contain ${BAND}"; return; fi
			done
		fi
		
		if [ "${CATALOG}" == "des_dr2" ]; then
			DES_BANDS=('g' 'r' 'i' 'z' 'y')
			DES_MAPS=('wavg_mag_psf_g' 'wavg_mag_psf_r' 'wavg_mag_psf_i' 'wavg_mag_psf_z' 'wavg_mag_psf_y')
			
			for i in ${!DES_BANDS[@]}; do
				if [ "${DES_BANDS[$i]}" == "$BAND" ]; then 
					jointcal_formatter "${CATALOG}" "${DES_MAPS[$i]}" "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: DES does not contain ${BAND}"; return; fi
			done
		fi

		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/jointcal_template.sh > ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s/photom_ref/${CATALOG}/g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh

		echo "Submitting to slurm..."
		sbatch ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sleep 20m
		
		# Soren's idea, submit successive jobs so that they wait in the queue for the previous to finish!
		# jointcal is single-threaded, so we don't need to hand it all resources, but I'll leave this here in case we need it later
		#
		# if [ -z "$JOBID" ]; then
		#	JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh)
		#	echo "Submitted jointcal_${BAND} with ${JOBID}"
		#else
		#	JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh)
		#	echo "Submitted jointcal_${BAND} with ${JOBID}"
		#fi

	done

}

# STEP 11: final visit summary, this does a little bit of data management and implements data from jointcal to produce finalized per-visit catalogs

# TODO (later): test different selection algorithms for aperture-correction...
# what matters for now is that enough ccd's have the correction so that, when they are stacked, the correction will average-out


# config-formatting helper function

final_visit_summary_formatter () {
	
	## text-fill fields ##
	
	# band
	PROCESS_BAND="${1}"
	
	# here we can add a check for the band and tweak these parameters to specific defaults
	# if we want to or need to...
	
	# psf stuff first
	PIFF_ORDER="${2}"
	PIFF_STAMP="${3}"
	
	# then the SN for aperture correction
	AP_CORR_SN="${4}"
	
	# running the actual text-replacement
	sed "s/piff_order/${PIFF_ORDER}/g" "${AUTO_PIPELINE_DIR}/config_templates/final_visit_summary_config_template.py" > "${CLUSTER_DIR}/configs/final_visit_summary_config_${PROCESS_BAND}.py"
	sed -i "s/piff_stamp/${PIFF_STAMP}/g" "${CLUSTER_DIR}/configs/final_visit_summary_config_${PROCESS_BAND}.py"
	sed -i "s/ap_corr_sn/${AP_CORR_SN}/g" "${CLUSTER_DIR}/configs/final_visit_summary_config_${PROCESS_BAND}.py"
	
}

final_visit_summary () {
	
	echo "Running STEP 11: final_visit_summary"

	for BAND in "$@"; do
		
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/final_visit_summary_template.sh > ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh
		# pass the current lsst_pipeline into the template
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh
		# pass the band being processed into the template
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh
		
		# by default I'll pass the LSP defaults to this config
		final_visit_summary_formatter "${BAND}" 2 25 200

		echo "Submitting final_visit_summary_${BAND} slurm..."
		
		# Soren's idea, submit successive jobs so that they wait in the queue for the previous to finish.
		# This let's us better utilize resources, at the cost of some overhead ran in serial at the beginning/end of each job
		if [ -z "$JOBID" ]; then
			JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh)
			echo "Submitted final_visit_summary_${BAND} with ${JOBID}"
		else
			JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh)
			echo "Submitted final_visit_summary_${BAND} with ${JOBID}"
		fi
		
		# sbatch ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh
		# sleep 5m
		
	done

	prompt_wait
}

# STEP 12: perform the first half of drp#step3a, (warp,assemble,detect,deblend)

coadd_3a () {

	echo "Running STEP 12: coadd_3a"
	
	# coadd_3a tends to have a lot of memory-related failures
	# so I've added this step to create a custom config for running 
	echo "Creating new BPS config..."
	
	QUERY_LEN=$(< query_result_${CLUSTER_NAME}.csv wc -l)
	MP=$((QUERY_LEN*60*4/(1000*5)))
	TASKS=$((RAM/MP))
	if [[ "${TASKS}" -gt "${CORES}" ]]; then
		ALLOC=$((CORES/NODES))
	else
		ALLOC=$((TASKS/NODES))
	fi
	MEM_ALLOC=$((MP*ALLOC))
	echo "Allocating ${NODES} nodes, ${ALLOC} cores, and ${MEM_ALLOC} GB of RAM per node."
	bps_config_formatter ${NODES} "${ALLOC}" "${MEM_ALLOC}" "48:00:00" "_coadd3a"
	
	# the configs below are now specified in the DRP-LoVoCCS.yaml directly
	# but they're here in case we need to tweak them more carefully
	# cp ${AUTO_PIPELINE_DIR}/config_templates/coadd_3a_assemble_config_template.py "${CLUSTER_DIR}/configs/coadd_3a_assemble_config.py"
	# cp ${AUTO_PIPELINE_DIR}/config_templates/coadd_3a_makewarp_config_template.py "${CLUSTER_DIR}/configs/coadd_3a_makewarp_config.py"
	# cp ${AUTO_PIPELINE_DIR}/config_templates/coadd_3a_deblend_config_template.py "${CLUSTER_DIR}/configs/coadd_3a_deblend_config.py"
	# cp ${AUTO_PIPELINE_DIR}/config_templates/coadd_3a_detection_config_template.py "${CLUSTER_DIR}/configs/coadd_3a_detection_config.py"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3a_template.sh > ${PROCESSING_STEP_DIR}/coadd_3a.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3a.sh
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/coadd_3a.sh
	prompt_wait

}

# STEP 13: coadd_3b, this step measures sources on the coadd

# config-formatting helper function

coadd_3b_formatter () {
	
	## text-fill fields ##
	
	# gaia is hard-coded as the default for astrometry
	ASTROM_REF="gaia"
	ASTROM_BAND="phot_g_mean_mag"
	
	# photometry is passed to the formatter
	PHOTOM_REF="${1}"
	PHOTOM_MAG="${2}"
	
	# last argument is the physical band being processed
	PROCESS_BAND="${3}"
	
	# running the actual text-replacement
	sed "s/astrom_ref/${ASTROM_REF}/g" ${AUTO_PIPELINE_DIR}/config_templates/coadd_3b_config_template.py > "${CLUSTER_DIR}/configs/coadd_3b_config_${PROCESS_BAND}.py"
	sed -i "s/astrom_band/${ASTROM_BAND}/g" "${CLUSTER_DIR}/configs/coadd_3b_config_${PROCESS_BAND}.py"
	sed -i "s/photom_ref/${PHOTOM_REF}/g" "${CLUSTER_DIR}/configs/coadd_3b_config_${PROCESS_BAND}.py"
	sed -i "s/photom_mag/${PHOTOM_MAG}/g" "${CLUSTER_DIR}/configs/coadd_3b_config_${PROCESS_BAND}.py"
	sed -i "s/process_band/${PROCESS_BAND}/g" "${CLUSTER_DIR}/configs/coadd_3b_config_${PROCESS_BAND}.py"
}

coadd_3b () {

	echo "Running STEP 13: coadd_3b"
	
	# coadd_3b used to be quite painful, but with the catalog-matching disabled it runs as one BIG step :)
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3b_template.sh > ${PROCESSING_STEP_DIR}/coadd_3b.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3b.sh
	
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/coadd_3b.sh
	prompt_wait

}

#TODO this step is now deprecated, remove after verifying that the new coadd_3a/coadd_3b/coadd_3c all work correctly
# STEP 14: finish-up drp#step3 (merge,forcedPhoto,consolidateTable)

coadd_3c () {

	echo "Running STEP 14: coadd_3c"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3c_template.sh > ${PROCESSING_STEP_DIR}/coadd_3c.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3c.sh
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/coadd_3c.sh
	prompt_wait

}

# STEP 15: produce a separate skycorr stack and detect sources on it

coadd_3c () {

	echo "Running STEP 15: coadd_3dc"

	# copying the config templates to the cluster-config folder, for now we have no custom configs to pass to these
	cp ${AUTO_PIPELINE_DIR}/config_templates/coadd_3c_skycorr_config_template.py "${CLUSTER_DIR}/configs/coadd_3c_skycorr_config.py"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3c_template.sh > ${PROCESSING_STEP_DIR}/coadd_3c.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3c.sh
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/coadd_3c.sh
	prompt_wait

}

# STEP 16: export data from the LSP and draw fov

export_data () {

	echo "Running STEP 16: export_data"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/export_data_template.sh > ${PROCESSING_STEP_DIR}/export_data.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/export_data.sh

	# create a mock-up of the gen2 directory structure w/ calexps and cats
	mkdir read_catalog_all_output
	mkdir combine_patch_color_output

	sbatch ${PROCESSING_STEP_DIR}/export_data.sh
}


# STEP 17: run photometric correction

#from here-out these are copies of SF's scripts
#TODO refractor these scripts to fit better in the "processing_step" format and check for bugs

photometric_correction () {

    REFCAT_INSTRUMENT=$1
    REFCAT_DATA_RELEASE=$2
    
    prompt "Running STEP 17: photometric_correction"

    sed "s/cluster_name/${CLUSTER_NAME}/g; s/photom_ref/${REFCAT_INSTRUMENT}/g; s/photom_dr/${REFCAT_DATA_RELEASE}/g" ${TEMPLATE_DIR}/photometric_correction_template.sh > ${PROCESSING_STEP_DIR}/photometric_correction.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/photometric_correction.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/photometric_correction.sh
    prompt_wait
    
    sleep 5s

    
}

# STEP 18: measure photometric redshifts

photo_z () {
    
    prompt "Running STEP 18: photo_z"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/photo_z_template.sh > ${PROCESSING_STEP_DIR}/photo_z.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/photo_z.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/photo_z.sh
    prompt_wait
    
    sleep 5s

    
}

# STEP 19: shear calibration

shear_calibration () {
    
    prompt "Running STEP 19: shear_calibration"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/shear_calibration_template.sh > ${PROCESSING_STEP_DIR}/shear_calibration.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/shear_calibration.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/shear_calibration.sh
    prompt_wait
    
    sleep 5s
 
}


# STEP 20: construct the mass map

mass_map () {
    
    prompt "Running STEP 20: mass_map"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/mass_map_template.sh > ${PROCESSING_STEP_DIR}/mass_map.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/mass_map.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/mass_map.sh
    prompt_wait
    
    sleep 5s
 
}

# STEP 21: fit the mass

mass_fit () {
    
    prompt "Running STEP 21: mass_fit"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/mass_fit_template.sh > ${PROCESSING_STEP_DIR}/mass_fit.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/mass_fit.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/mass_fit.sh
    prompt_wait
    
    sleep 5s
 
}

# STEP 22: check quality

#TODO We can use the data saved from check/select_visit to make a much more accurate map of the depth and pointings across the patches
quality_check () {
    
    prompt "Running STEP 22: quality_check"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/quality_check_template.sh > ${PROCESSING_STEP_DIR}/quality_check.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/quality_check.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/quality_check.sh
    prompt_wait
    
    sleep 5s
 
}

# STEP 23: red sequence distribution

red_sequence () {
    
    prompt "Running STEP 23: red_sequence"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/red_sequence_template.sh > ${PROCESSING_STEP_DIR}/red_sequence.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/red_sequence.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/red_sequence.sh
    prompt_wait
    
    sleep 5s
 
}

# STEP 24: meta_4a

meta_4a () {

	echo "Running STEP 24: meta_4a"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta4a_template.sh > ${PROCESSING_STEP_DIR}/meta_4a.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_4a.sh
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/meta_4a.sh
	prompt_wait

}

# STEP 25: meta_4b

meta_4b () {

	echo "Running STEP 26: meta_4b"
	
	
	# copying over the config templates, for now these are customized
	cp ${AUTO_PIPELINE_DIR}/config_templates/meta_4b_forced_config_template.py "${CLUSTER_DIR}/configs/meta_4b_forced_config.py"
	cp ${AUTO_PIPELINE_DIR}/config_templates/meta_4b_measure_config_template.py "${CLUSTER_DIR}/configs/meta_4b_measure_config.py"
	
	for SHEARTYPE in "noshear" "1p" "1m" "2p" "2m"; do
	
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta4b_template.sh > ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh
		
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh
		sed -i "s/shear_type/${SHEARTYPE}/g" ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh
	
		echo "Submitting to slurm..."
	
		if [ -z "$JOBID" ]; then
			JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh)
			echo "Submitted meta_4b_${SHEARTYPE} with ${JOBID}"
		else
			JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh)
			echo "Submitted meta_4b_${SHEARTYPE} with ${JOBID}"
		fi

	done
	
	prompt_wait

}

# STEP 26: export metadetect data

meta_export () {

	echo "Running STEP 26: meta_export"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_export_template.sh > ${PROCESSING_STEP_DIR}/meta_export.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_export.sh

	# create a mock-up of the gen2 directory structure w/ calexps and cats
	mkdir metadetect_export

	sbatch ${PROCESSING_STEP_DIR}/meta_export.sh
}

# STEP 27: process metadetect data

meta_processing () {

	echo "Running STEP 27: meta_processing"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_processing_template.sh > ${PROCESSING_STEP_DIR}/meta_processing.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_processing.sh

	mkdir metadetect_processing

	sbatch ${PROCESSING_STEP_DIR}/meta_processing.sh
}

# STEP 28: lensing w. metadetect

meta_lensing () {

	echo "Running STEP 28: meta_processing"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_lensing_template.sh > ${PROCESSING_STEP_DIR}/meta_lensing.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_lensing.sh

	sbatch ${PROCESSING_STEP_DIR}/meta_lensing.sh
}

# STEP 29: gotta blast

gotta_blast () {

    xmin=$1
    ymin=$2
    xmax=$3
    ymax=$4
    
    prompt "Running STEP 29: gotta_blast"
    
    read -p "Are you sure? Deletion cannot be undone!" -n 1 -r
    echo    # (optional) move to a new line
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
    fi
    
    sed "s/cluster_name/${CLUSTER_NAME}/g; s/xmin/${xmin}/g; s/ymin/${ymin}/g; s/xmax/${xmax}/g; s/ymax/${ymax}/g" ${TEMPLATE_DIR}/gotta_blast_template.sh > ${PROCESSING_STEP_DIR}/gotta_blast.sh

    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/gotta_blast.sh

    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/gotta_blast.sh
    prompt_wait
    
    sleep 5s
 
}

# == RUNNING STEPS == #

#STEP0 NOTES: 
# After copying and pasting run_steps_Gen3 into a directory with the Cluster Name, create_output creates a series of folders and python scripts. This only takes a few seconds to run

create_output


#STEP1 NOTES: 
# This step queries and downloads frames from the NoirLab Science Archive
# it usually takes a few hours too run, but if there is a lot of traffic or
# other bandwidth issues (see the TODO on line 68) it can take as long as 12 hrs

#download_raw


#STEP2 NOTES:
# Unfortunately, LSSTPipe has doesn't have any code built to catch errors due to
# missing bits in the fits header or other issues which may corrupt a fits file.
# As a result we have a step which manually opens checks each file (by opening each file python
# and checking to see if the script crashes dramatically). This should only take a few minutes to run

#check_raws


#STEP3 NOTES: 
# Really... STEP2 is one script which I wrapped into a job-array to help things run a bit faster
# but because of how it was written, the corrupt fits have to be moved to a separate directory in another script.
# Eventually I'll wrap these into one step..... but I don't want to manually rename all the step numbers below
# at the moment.

#move_corrupt_raws


#STEP4 NOTES:
# To begin LSSTPipe, we need to create a repository which will manage all of the files produced during processing
# this is a quick process, so like the above steps this is run from the command line rather than with a batch script

#initialize_repo


#STEP5 NOTES:
# Next, we have to ingest our data into the repository. This step takes care of ingesting all of the catalogs,
# calibrations, and raw data into the repo. Depending on the amount of raw data, this can take an hour at most
# but usually takes less.

#ingest_data


#STEP6 NOTES: 
# This technically runs step1 of the Data Release Pipe (DRP), the big steps this includes are:
# ISR (Instrumental Signal Removal: Crosstalk, Nonlinearity, Bias, Fringe, Flat, Brighter-Fatter)
# Image Characterization (Initial Measurements, Cosmic-Ray Repair, Background Subtraction, PSF Measured)
# Image Calibration (Initial Astrometry & Photometry)
# Remaining tasks consolidate sources into per-detector tables

# For photometric refcats, the options are des_dr2 (g,r,i,z,y), ps1 (g,r,i,z,y), sdss (u-band only), sm_dr1_(g,i,r,v,z), sm_dr2_BAND_(g,i,r,v,z), see ~/Clusters/calib_catalog_repo/catalogs_new/CLUSTER for source-catalog coverage.
# This is the first BIG processing step, which takes a VERY long time to run, typically on the order of ~12hrs
# (per-band, 20 cores allocated) so try to run it overnight if you can.

# Occassionally, these separate jobs will happen to write to the repo at the same time, which will cause a failure
# When this occurs, just rerun the corresponding band using sbatch...

# Also very rarely when a job is submitted although the bps.USRN job will be created,
# it will not be writing any outputs. In this case the fix is simply cancelling that job...
# The "admin" script (process_ccd_BAND_CLN) will automatically resubmit a new set of managers/workers.

#process_ccd sdss,u ps1,g ps1,r ps1,i ps1,z ps1,y


#STEP7+8 NOTES: 
# Before moving forward, we check the seeing in each visit/ccd and trim visits which are not "lensing-quality"
# We fit a Moffat profile to reference stars for measuring the FWHM and use the second moments to measure the distortion
# in the psf. Technically there is a function for doing this embedded in LSSTPipe, but we have yet to test it. After selecting
# the good detectors, a skymap object is created which encloses them
# These usually take a few minutes to run each...

#check_visit u g i r z
#select_visit


#STEP9 NOTES:
# This runs step2a of DRP. Taking the best CCD's following select_visit, we create final visit summaries
# which are required for jointcal. This only takes a few minutes to run...

#visit_summary u g r i z


#STEP10 NOTES: 
# This runs step2b of DRP, joint-calibration. For precision astrophysics, errors in photometry/astrometry must be minimized.
# It happens that, on average, you can achive a more precise calibration by fitting the astrometry/photometry of each visit
# individually, rather than matching to a reference catalog after stacking. This can be done by fully modelling the
# atmosphere/optics (fgcm) or using an emiprical correction which models variations in brightnesses and positions
# between exposures (jointcal). Jointcal models variations in positions/brightnesses with a polynomial anchored in
# refcats and seeks to minimize a joint chi-squared. Ususally this takes about an hour/band (20-cores)...

# Options are des_dr2 (g,r,i,z,y), ps1 (g,r,i,z,y), sdss (u-band only), sm_dr1_(g,i,r,v,z), sm_dr2_BAND_(g,i,r,v,z) the arguments here and in coadd_3b should be identical to the input of process_ccd

#jointcal sdss,u ps1,g ps1,r ps1,i ps1,z


#STEP11 NOTES: 
# This runs step2d of DRP. For reasons I can't quite remember, step2c is optional and isn't currently functional on DECam.
# It includes a final round of calibration which applies the corrections derived from jointcal and creates finalized
# visit summaries. This takes ~1hr/band with 20-cores

#final_visit_summary u g r i z


#STEP12-15 NOTES: 
# These steps run through step3a,step3b,step3c, and step3d of DRP.
# step3a runs coaddition, detects sources, and runs the deblender (~3 hours w/ 140-cores)
# step3b runs measurements on the coadds (~4 hrs/band w/ 20-cores)
# step3c consolidates everything into tables and runs forced photometry (~5-7 hours w/ 140-cores)
# step3d runs skycorrection, then makes a skycorrected coadded stack (~1 hour w/ 140-cores)
# Technically 3a,3b,3c can be run together if we disable the refcat matching... but since 3b is prone to errors
# it may not be very useful to merge them

#coadd_3a
#coadd_3b sdss,u ps1,g ps1,r ps1,i ps1,z
#coadd_3c
#coadd_3d

#STEP16 NOTES:
# This step exports all of the data out of LSSTPipe, including:
# Object catalogs containing shapes and magnitudes
# Fits images containing the full fov
# irg-images displaying the data
# roughly ~1hr to run

#export_data


# == Processing with LSSTPipe is mostly done and scripts only take a few minutes to run from here on == #


#STEP17 NOTES:
# From this point forward, the analysis is carried out on a catalog-level, the first step to this is
# correcting the photometry for color-terms and extinction (called de-redding). u-band is particularly funky here since
# there are very few existing refcats to calibrate the data with... so instead we use mock-observations
# created by integrating the spectra of reference-stars with the known transmission of the DECam filters
# and measure the u-band color terms directly from that calculation ("stellar-locus correction").

# Use "ps1" "dr1" or "sm" "dr2" depending on the refcat used during process_ccd.

#photometric_correction "ps1" "dr1"


#STEP18 NOTES:
# Once we have a corrected catalog, it's time to identify galaxies in the field and estimate their redshift.
# For this we use a Bayesian PhotoZ algorithm.

#photo_z


#STEP19 NOTES:
# For the time being, we run a naive shear-calibration (e.g. g ~ e/(2R) without any weights).
# This is almost certainly not the best we can do, but based on Soren's last Source-Injection
# batch it's certainly better than usign the HSC method

#shear_calibration


#STEP20 NOTES:
# Now that we have photo-z's, we can transform the HSM shape measurements into reduced shears and assemble a 
# mass_map. Generally, what we produce at this stage is not directly contours representing the mass, but instead
# a "mass-aperture map", which convolves the shear field with a filter built to pick out a particular signal.
# Our filter is built to be particularly sensitive to NFW-like structures, which will have a high SN-ratio

#mass_map


#STEP21 NOTES:
# This step actually extracts the mass by fitting the shears of an NFW-profile to the observed shear-field and
# bootstrapping to produce error-bars.

#mass_fit


#STEP22 NOTES:
# This stage looks for correlations between different objects and produces plots summarizing the quality of 
# observations/data which was used.

#quality_check


#STEP23 NOTES:
# Generally, in a cluster the "red-sequence galaxies" are the oldest objects in the cluster and, as a result,
# effectively trace out the distribution of mass in the cluster. The red-sequence is selected from an HR-diagram
# of the cluster and smoothed contours representing their number density are produced

#red_sequence


# == Briefly going back to LSSTPipe for running metadetect! == #


#STEP24-28 NOTES:
# these steps run our implementation of metadetect
# meta_4a runs the shears on our coadds
# and meta_4b runs detect/deblend/measure on them, which we can use to build a robust calibration
# meta_export collects all of the outputs from the 5 different versions of our coadds
# meta_processing/meta_lensing finally process those coadds and run the lensing portion (including shear-calibration!)

#meta_4a
#meta_4b
#meta_export
#meta_processing
#meta_lensing


#STEP29 NOTES:
# blast intermediate collections and other datasets, then zip the submit-directory

#gotta_blast


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

    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/triaxiality.sh

    echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/triaxiality.sh
    prompt_wait
    
    sleep 5s
 
}








