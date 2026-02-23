# Processing pipeline for LoVoCCS!

# ==== #

# This script is built to work users through each task for processing clusters for the Local Volume Complete Cluster Survey (LoVoCCS)
# See the readme for more information!

# == Name of the cluster to process == #
# this may need to be changed for a few clusters...
CLUSTER_NAME='A85'


# == Resources to allocate for cluster processing == #
# Check out https://ccv.brown.edu/rates/ to see how many cores and how much RAM you have available 
# If you're not on Oscar, hopefully you already know what this looks like
# Make sure to leave some cores/memory available so you can load up a terminal or do other science
NODES=10 
CORES=180 # multiple of NODES only!
RAM=1300 # multiple of NODES only!
WALL_TIME=48:00:00 


# == VARIABLES == #
# This is setup for folks operating on our Oscar/CCV allocation
# Users outside of Brown will need to tweak this
AUTO_PIPELINE_DIR="/gpfs/data/idellant/Clusters/gen3_processing/lovoccs_pipe" # path to your installation of lovoccs_pipe
TEMPLATE_DIR="${AUTO_PIPELINE_DIR}/processing_step_templates" # shouldn't need to be adjusted
LOAD_PIPELINE_PATH="/gpfs/data/idellant/lsst_stack_v26_0_0/loadLSST.bash" # path to your v26.0.0 LSP installation
CLUSTER_DIR="/gpfs/data/idellant/Clusters/gen3_processing/${CLUSTER_NAME}" # the directory in which you want to process
PROCESSING_STEP_DIR="${CLUSTER_DIR}/processing_step" # shouldn't need to be adjusted
CALIB_CATALOG_REPO="/gpfs/data/idellant/Clusters/calib_catalog_repo" # location of calibrations and refcats


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

download_raw_format () {

	echo "Running STEP 1: download_raw"

	echo "Downloading ${CLUSTER_NAME}..."
	echo "Raw images will be downloaded into .../${CLUSTER_NAME}/raw"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/noao_download_manager_template.sh > ${PROCESSING_STEP_DIR}/noao_download_manager.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/noao_download_manager.sh

	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/noao_download_manager.sh
	#prompt_wait

}


# STEP 2: Check for any exposures with a bad/corrupt fits header

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

check_raws_format () {

	echo "Running STEP 2: check_raws"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/check_raws_template.sh > ${PROCESSING_STEP_DIR}/check_raws.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/check_raws.sh

	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/check_raws.sh
	#prompt_wait

}


# STEP 3: Move those exposures to a separate folder

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

move_corrupt_raws_format () {

	echo "Running STEP 3: move_corrupt_raws"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/move_corrupt_raws_template.sh > ${PROCESSING_STEP_DIR}/move_corrupt_raws.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/move_corrupt_raws.sh

	#echo "Now running..."
	#bash ${PROCESSING_STEP_DIR}/move_corrupt_raws.sh
	#prompt_wait

}


# STEP 4: initialize the repository & create the Butler

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
	butler write-curated-calibrations repo/repo lsst.obs.decam.DarkEnergyCamera

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

ingest_data_format () {

	echo "Running STEP 5: ingest_data"

	# pass the cluster name into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/ingest_data_template.sh > ${PROCESSING_STEP_DIR}/ingest_data.sh

	# pass the current lsst_pipeline and cluster_dir into the script
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/ingest_data.sh

	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/ingest_data.sh
	#prompt_wait

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

# config-formatting helper function
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

# the actual processing-step
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
		
		# sky-mapper dr4 formatting
		if [ "${CATALOG}" == "sm_dr4" ]; then
			SM4_BANDS=('u' 'g' 'r' 'i' 'z')
			SM4_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM4_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM4_BANDS[@]}; do
				if [ "${SM4_BANDS[$i]}" == "$BAND" ]; then 
					process_ccd_refcat_formatter "${CATALOG}_${SM4_CAT_BANDS[$i]}" "${SM4_MAPS[$i]}" "${BAND}"
					process_ccd_psf_formatter "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR4 does not contain ${BAND}"; return; fi
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

		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/process_ccd_template.sh > ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh

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

process_ccd_format () {

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
		
		# sky-mapper dr4 formatting
		if [ "${CATALOG}" == "sm_dr4" ]; then
			SM4_BANDS=('u' 'g' 'r' 'i' 'z')
			SM4_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM4_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM4_BANDS[@]}; do
				if [ "${SM4_BANDS[$i]}" == "$BAND" ]; then 
					process_ccd_refcat_formatter "${CATALOG}_${SM4_CAT_BANDS[$i]}" "${SM4_MAPS[$i]}" "${BAND}"
					process_ccd_psf_formatter "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR4 does not contain ${BAND}"; return; fi
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

		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/process_ccd_template.sh > ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh

		# Soren's idea, submit successive jobs so that they wait in the queue for the previous to finish.
		# This let's us better utilize resources, at the cost of some overhead ran in serial at the beginning/end of each job
		#if [ -z "$JOBID" ]; then
		#	JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh)
		#	echo "Submitted process_ccd_${BAND} with ${JOBID}"
		#else
		#	JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/process_ccd_${BAND}.sh)
		#	echo "Submitted process_ccd_${BAND} with ${JOBID}"
		#fi
	done

}

# STEP 7: check each visit to measure the seeing (PSF-FWHM) and average ellipticity across the sensor

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

check_visit_format () {
	
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
		#sbatch ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		
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
		#sbatch ${PROCESSING_STEP_DIR}/check_visit_${BAND}.sh
		
		;;
		*)

		echo "${BAND} is not allowed... please choose between ugirzy!"
		;;
		esac
	done

	#prompt_wait
}


# STEP 8: apply star-count/fwhm/ellip cuts and create a skymap enclosing our exposures

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

select_visit_format () {

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

	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/select_visit.sh
	#prompt_wait

}

# STEP 9: consolidate tables from process_ccd and ready data for jointcal

visit_summary () {

	echo "Running STEP 9: visit_summary"
	
	for BAND in "$@"; do
		
		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/visit_summary_template.sh > ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		
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

visit_summary_format () {

	echo "Running STEP 9: visit_summary"
	
	for BAND in "$@"; do
		
		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/visit_summary_template.sh > ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh
		
		# Soren's idea, submit successive jobs so that they wait in the queue for the previous to finish.
		# This let's us better utilize resources, at the cost of some overhead ran in serial at the beginning/end of each job
		#if [ -z "$JOBID" ]; then
		#	JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh)
		#	echo "Submitted visit_summary_${BAND} with ${JOBID}"
		#else
		#	JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/visit_summary_${BAND}.sh)
		#	echo "Submitted visit_summary_${BAND} with ${JOBID}"
		#fi


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
		JOINTCAL_ORDER=2 # order 2 is sufficient!
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
	
	# delay successive submissions by ~10min to prevent conflicts later
	# this isn't foolproof... but works surprisingly well!
	DELAY=0
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
		
		# sky-mapper dr4 formatting
		if [ "${CATALOG}" == "sm_dr4" ]; then
			SM4_BANDS=('u' 'g' 'r' 'i' 'z')
			SM4_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM4_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM4_BANDS[@]}; do
				if [ "${SM4_BANDS[$i]}" == "$BAND" ]; then 
					jointcal_formatter "${CATALOG}_${SM4_CAT_BANDS[$i]}" "${SM4_MAPS[$i]}" "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR4 does not contain ${BAND}"; return; fi
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

		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/jointcal_template.sh > ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s/photom_ref/${CATALOG}/g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh

		echo "Submitting to slurm..."
		
		sbatch --begin=now+${DELAY} ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		DELAY=$((${DELAY} + 600))

	done

}

jointcal_format () {

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
	
	# delay successive submissions by ~10min to prevent conflicts later
	# this isn't foolproof... but works surprisingly well!
	DELAY=0
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
		
		# sky-mapper dr4 formatting
		if [ "${CATALOG}" == "sm_dr4" ]; then
			SM4_BANDS=('u' 'g' 'r' 'i' 'z')
			SM4_MAPS=('v_psf' 'g_psf' 'r_psf' 'i_psf' 'z_psf')
			SM4_CAT_BANDS=('v' 'g' 'r' 'i' 'z')
			
			for i in ${!SM4_BANDS[@]}; do
				if [ "${SM4_BANDS[$i]}" == "$BAND" ]; then 
					jointcal_formatter "${CATALOG}_${SM4_CAT_BANDS[$i]}" "${SM4_MAPS[$i]}" "${BAND}"
					break
				fi
				if [ $i == 4 ]; then echo "ERROR: SM_DR4 does not contain ${BAND}"; return; fi
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

		# pass the cluster name and band into the template
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/jointcal_template.sh > ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s/process_band/${BAND}/g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s/photom_ref/${CATALOG}/g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh

		#echo "Submitting to slurm..."
		
		#sbatch --begin=now+${DELAY} ${PROCESSING_STEP_DIR}/jointcal_${BAND}.sh
		#DELAY=$((${DELAY} + 600))

	done

}


# STEP 11: final visit summary, this does a little bit of data management and implements data from jointcal to produce finalized per-visit catalogs

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
		
	done

	prompt_wait
}

final_visit_summary_format () {
	
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
		#if [ -z "$JOBID" ]; then
		#	JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh)
		#	echo "Submitted final_visit_summary_${BAND} with ${JOBID}"
		#else
		#	JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/final_visit_summary_${BAND}.sh)
		#	echo "Submitted final_visit_summary_${BAND} with ${JOBID}"
		#fi
		
	done

	prompt_wait
}


# STEP 12: perform the first half of drp#step3a, (warp,assemble)
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

coadd_3a_format () {

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
	
	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/coadd_3a.sh
	#prompt_wait

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

	# pass the cluster name and band into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3b_template.sh > ${PROCESSING_STEP_DIR}/coadd_3b.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3b.sh
	echo "Submitting to slurm..."
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/coadd_3b.sh
	prompt_wait

}

coadd_3b_format () {

	echo "Running STEP 13: coadd_3b"

	# pass the cluster name and band into the template
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3b_template.sh > ${PROCESSING_STEP_DIR}/coadd_3b.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3b.sh
	echo "Submitting to slurm..."
	
	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/coadd_3b.sh
	#prompt_wait

}


# STEP 14: finish-up drp#step3 (merge,forcedPhoto,consolidateTable)
coadd_3c () {

	echo "Running STEP 14: coadd_3c"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3c_template.sh > ${PROCESSING_STEP_DIR}/coadd_3c.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3c.sh
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/coadd_3c.sh
	prompt_wait

}

coadd_3c_format () {

	echo "Running STEP 14: coadd_3c"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3c_template.sh > ${PROCESSING_STEP_DIR}/coadd_3c.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3c.sh
	
	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/coadd_3c.sh
	#prompt_wait

}


# STEP 15: produce a separate skycorr stack and detect sources on it
coadd_3d () {

	echo "Running STEP 15: coadd_3d"

	# copying the config templates to the cluster-config folder, for now we have no custom configs to pass to these
	cp ${AUTO_PIPELINE_DIR}/config_templates/coadd_3d_skycorr_config_template.py "${CLUSTER_DIR}/configs/coadd_3d_skycorr_config.py"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3d_template.sh > ${PROCESSING_STEP_DIR}/coadd_3d.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3d.sh
	
	echo "Submitting to slurm..."
	sbatch ${PROCESSING_STEP_DIR}/coadd_3d.sh
	prompt_wait

}

coadd_3d_format () {

	echo "Running STEP 15: coadd_3d"

	# copying the config templates to the cluster-config folder, for now we have no custom configs to pass to these
	cp ${AUTO_PIPELINE_DIR}/config_templates/coadd_3d_skycorr_config_template.py "${CLUSTER_DIR}/configs/coadd_3d_skycorr_config.py"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/coadd3d_template.sh > ${PROCESSING_STEP_DIR}/coadd_3d.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/coadd_3d.sh
	
	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/coadd_3d.sh
	#prompt_wait

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

export_data_format () {

	echo "Running STEP 16: export_data"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/export_data_template.sh > ${PROCESSING_STEP_DIR}/export_data.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/export_data.sh

	# create a mock-up of the gen2 directory structure w/ calexps and cats
	mkdir read_catalog_all_output
	mkdir combine_patch_color_output

	#sbatch ${PROCESSING_STEP_DIR}/export_data.sh
}


# STEP 17: run photometric correction
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

photometric_correction_format () {

    REFCAT_INSTRUMENT=$1
    REFCAT_DATA_RELEASE=$2
    
    prompt "Running STEP 17: photometric_correction"

    sed "s/cluster_name/${CLUSTER_NAME}/g; s/photom_ref/${REFCAT_INSTRUMENT}/g; s/photom_dr/${REFCAT_DATA_RELEASE}/g" ${TEMPLATE_DIR}/photometric_correction_template.sh > ${PROCESSING_STEP_DIR}/photometric_correction.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/photometric_correction.sh
    
    #echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/photometric_correction.sh
    #prompt_wait
    
    #sleep 5s
 
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

photo_z_format () {
    
    prompt "Running STEP 18: photo_z"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/photo_z_template.sh > ${PROCESSING_STEP_DIR}/photo_z.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/photo_z.sh
    
    #echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/photo_z.sh
    #prompt_wait
    
    #sleep 5s

    
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

shear_calibration_format () {
    
    prompt "Running STEP 19: shear_calibration"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/shear_calibration_template.sh > ${PROCESSING_STEP_DIR}/shear_calibration.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/shear_calibration.sh
    
    #echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/shear_calibration.sh
    #prompt_wait
    
    #sleep 5s
 
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

mass_map_format () {
    
    prompt "Running STEP 20: mass_map"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/mass_map_template.sh > ${PROCESSING_STEP_DIR}/mass_map.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/mass_map.sh
    
    #echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/mass_map.sh
    #prompt_wait
    
    #sleep 5s
 
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

mass_fit_format () {
    
    prompt "Running STEP 21: mass_fit"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/mass_fit_template.sh > ${PROCESSING_STEP_DIR}/mass_fit.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/mass_fit.sh
    
    #echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/mass_fit.sh
    #prompt_wait
    
    #sleep 5s
 
}


# STEP 22: check quality
quality_check () {
    
    prompt "Running STEP 22: quality_check"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/quality_check_template.sh > ${PROCESSING_STEP_DIR}/quality_check.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/quality_check.sh
    
    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/quality_check.sh
    prompt_wait
    
    sleep 5s
 
}

quality_check_format () {
    
    prompt "Running STEP 22: quality_check"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/quality_check_template.sh > ${PROCESSING_STEP_DIR}/quality_check.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/quality_check.sh
    
    #echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/quality_check.sh
    #prompt_wait
    
    #sleep 5s
 
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

red_sequence_format () {
    
    prompt "Running STEP 23: red_sequence"

    sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/red_sequence_template.sh > ${PROCESSING_STEP_DIR}/red_sequence.sh
    
    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/red_sequence.sh
    
    #echo "Submitting sbatch script..."
    #sbatch ${PROCESSING_STEP_DIR}/red_sequence.sh
    #prompt_wait
    
    #sleep 5s
 
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

meta_4a_format () {

	echo "Running STEP 24: meta_4a"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta4a_template.sh > ${PROCESSING_STEP_DIR}/meta_4a.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_4a.sh
	
	#echo "Submitting to slurm..."
	#sbatch ${PROCESSING_STEP_DIR}/meta_4a.sh
	#prompt_wait

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

meta_4b_format () {

	echo "Running STEP 26: meta_4b"
	
	
	# copying over the config templates, for now these are customized
	cp ${AUTO_PIPELINE_DIR}/config_templates/meta_4b_forced_config_template.py "${CLUSTER_DIR}/configs/meta_4b_forced_config.py"
	cp ${AUTO_PIPELINE_DIR}/config_templates/meta_4b_measure_config_template.py "${CLUSTER_DIR}/configs/meta_4b_measure_config.py"
	
	for SHEARTYPE in "noshear" "1p" "1m" "2p" "2m"; do
	
		sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta4b_template.sh > ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh
		
		sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh
		sed -i "s/shear_type/${SHEARTYPE}/g" ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh
	
		#echo "Submitting to slurm..."
	
		#if [ -z "$JOBID" ]; then
		#	JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh)
		#	echo "Submitted meta_4b_${SHEARTYPE} with ${JOBID}"
		#else
		#	JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/meta_4b_${SHEARTYPE}.sh)
		#	echo "Submitted meta_4b_${SHEARTYPE} with ${JOBID}"
		#fi

	done
	
	#prompt_wait

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

meta_export_format () {

	echo "Running STEP 26: meta_export"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_export_template.sh > ${PROCESSING_STEP_DIR}/meta_export.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_export.sh

	# create a mock-up of the gen2 directory structure w/ calexps and cats
	mkdir metadetect_export

	#sbatch ${PROCESSING_STEP_DIR}/meta_export.sh
}


# STEP 27: process metadetect data
meta_processing () {

	echo "Running STEP 27: meta_processing"
	
	mkdir metadetect_processing
	
	# split this into 5 separate shears to process
	for SHEARTYPE in "noshear" "1p" "2p" "1m" "2m"; do 
	
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_processing_template.sh > ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh
	sed -i "s|shear_type|${SHEARTYPE}|g" ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh

	sbatch ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh
	
	done
	
	prompt_wait
	
}

meta_processing_format () {

	echo "Running STEP 27: meta_processing"
	
	mkdir metadetect_processing
	
	# split this into 5 separate shears to process
	for SHEARTYPE in "noshear" "1p" "2p" "1m" "2m"; do 
	
	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_processing_template.sh > ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh
	sed -i "s|shear_type|${SHEARTYPE}|g" ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh

	#sbatch ${PROCESSING_STEP_DIR}/meta_processing_${SHEARTYPE}.sh
	
	done
	
	#prompt_wait
	
}


# STEP 28: lensing w. metadetect
meta_lensing () {

	echo "Running STEP 28: meta_lensing"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_lensing_template.sh > ${PROCESSING_STEP_DIR}/meta_lensing.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_lensing.sh

	sbatch ${PROCESSING_STEP_DIR}/meta_lensing.sh
}

meta_lensing_format () {

	echo "Running STEP 28: meta_lensing"

	sed "s/cluster_name/${CLUSTER_NAME}/g" ${TEMPLATE_DIR}/meta_lensing_template.sh > ${PROCESSING_STEP_DIR}/meta_lensing.sh
	sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g;s|py_scripts|${AUTO_PIPELINE_DIR}/python_scripts|g" ${PROCESSING_STEP_DIR}/meta_lensing.sh

	#sbatch ${PROCESSING_STEP_DIR}/meta_lensing.sh
}



# STEP 29: gotta blast
gotta_blast () {

    xmin=$1
    ymin=$2
    xmax=$3
    ymax=$4
    
    prompt "Running STEP 24: gotta_blast"
    
    read -p "Are you sure? Deletion cannot be undone!" -n 1 -r
    echo    # (optional) move to a new line
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
    fi
    
    sed "s/cluster_name/${CLUSTER_NAME}/g; s/xmin/${xmin}/g; s/ymin/${ymin}/g; s/xmax/${xmax}/g; s/ymax/${ymax}/g" ${TEMPLATE_DIR}/gotta_blast_template.sh > ${PROCESSING_STEP_DIR}/gotta_blast.sh

    sed -i "s|load_pipeline_path|${LOAD_PIPELINE_PATH}|g;s|cluster_dir|${CLUSTER_DIR}|g" ${PROCESSING_STEP_DIR}/gotta_blast.sh

    echo "Submitting sbatch script..."
    sbatch ${PROCESSING_STEP_DIR}/gotta_blast.sh
    prompt_wait
    
    sleep 5s
 
}
