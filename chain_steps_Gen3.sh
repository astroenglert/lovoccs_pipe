# A modified version of run_steps which submits steps that typically succeed in rapid succession as a "chain"
# that sits in the queue and patiently waits for the previous job to succeed
# first source Gen3 functions
source functions_Gen3.sh

# CHAIN 1
# runs the download and removes any corrupt .fits files
chain_1 () {

# format all the bash scripts
download_raw_format
check_raws_format
move_corrupt_raws_format

# now submit them in a dependency chain
for TASK in "download_raw" "check_raws" "move_corrupt_raws"; do
	if [ -z "$JOBID" ]; then
		JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	else
		JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	fi
done

}


# ingest, process_ccd, check_visit
# as jointcal can fail occassionally
chain_2 () {

# format bash scripts
ingest_data_format
# EDITME
process_ccd_format sdss,u ps1,g ps1,r ps1,i ps1,z ps1,y
# EDITME
check_visit_format u g r i z

# now submit with the proper dependencies
for TASK in "ingest_data" "process_ccd_u" "process_ccd_g"  "process_ccd_r"  "process_ccd_i"  "process_ccd_z"  "process_ccd_y"; do
	if [ -z "$JOBID" ]; then
		JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	else
		JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	fi
done

# these can neatly run in parallel 
for TASK in "check_visit_u" "check_visit_g" "check_visit_r" "check_visit_i" "check_visit_z"; do
	JOBIDNEW=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/${TASK}.sh)
	echo "Submitted ${TASK} with ${JOBIDNEW}"
done

}

# once check_visit has its outputs validated, select_visit prunes for quality and visit_summary gets stuff ready for jointcal
chain_3 () {

select_visit_format
# EDITME
visit_summary_format u g r i z

for TASK in "select_visit" "visit_summary_u" "visit_summary_g" "visit_summary_r" "visit_summary_i" "visit_summary_z"; do
	if [ -z "$JOBID" ]; then
		JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	else
		JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	fi
done

}

# the coadd steps tend to have failures due to memory allocations
# so those are not chained together
# but all the post-LSP steps can be chained together
chain_4 () {

export_data_format
# EDITME
photometric_correction_format "ps1" "dr2"
photo_z_format
shear_calibration_format
mass_map_format
mass_fit_format
quality_check_format
red_sequence_format

for TASK in "export_data" "photometric_correction" "photo_z" "shear_calibration" "mass_map" "mass_fit" "quality_check" "red_sequence"; do
	if [ -z "$JOBID" ]; then
		JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	else
		JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	fi
done

}

# meta_4b tends to have failures as well, so we avoid chaining meta_4a/meta_4b
# but the rest of the meta steps can be safely chained
chain_5 () {

meta_export_format
meta_processing_format
meta_lensing_format

for TASK in "meta_export" "meta_processing" "meta_lensing"; do
	if [ -z "$JOBID" ]; then
		JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	else
		JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	fi
done

}

# if you want to take a chance and submit the coadd steps in a chain, you can
# but I won't list it below as generally these steps have OOM errors that need some careful management
chain_coadd () {

# format all the bash scripts
coadd_3a_format
coadd_3b_format
coadd_3c_format
coadd_3d_format

# now submit them in a dependency chain
for TASK in "coadd_3a" "coadd_3b" "coadd_3c" "coadd_3d"; do
	if [ -z "$JOBID" ]; then
		JOBID=$(sbatch --parsable ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	else
		JOBID=$(sbatch --parsable --dependency=afterany:${JOBID} ${PROCESSING_STEP_DIR}/${TASK}.sh)
		echo "Submitted ${TASK} with ${JOBID}"
	fi
done

}

# when using chains, the refcats/bands have to be edited above
# line numbers are provided to point you to the code that you need to edit
# otherwise its the same as run-steps, comment and un-comment as needed to run each step
# some of them have now been chained together to help reduce our idle time

# == STEP0 == #
#create_output


# == CHAIN 1 NOTES == #
# download_raw, check_raws, move_corrupt_raws
# LINES 8-26

# == CHAIN1 COMMAND == #
#chain_1


# == STEP4 == #
#initialize_repo


# == CHAIN2 NOTES == #
# ingest_data, process_ccd, check_visit
# LINES 35-38

# == CHAIN2 COMMAND == #
#chain_2


# == CHAIN3 NOTES == #
# select_visit, visit_summary
# LINES 62-63
# make sure to inspect check_visit/plots to ensure nothing has gone horribly wrong

# == CHAIN3 COMMAND == #
#chain_3


# == STEP10 == #
#jointcal sm_dr4,u ps1,g ps1,r ps1,i ps1,z


# == STEP11 == #
#final_visit_summary u g r i z


# == STEP12-15 == #
#coadd_3a
#coadd_3b
#coadd_3c
#coadd_3d


# == CHAIN4 NOTES == #
# export_data, photometric_correction, photo_z, shear_calibration, mass_fit, mass_map, red_sequence, quality_check
# LINES 83-84

# == CHAIN4 COMMAND == #
#chain_4


# == STEP24-25 == #
#meta_4a
#meta_4b


# == CHAIN5 NOTES == #
# meta_export, meta_processing, meta_lensing

# == CHAIN5 COMMAND == #
#chain_5


# == STEP29 == #
#gotta_blast









