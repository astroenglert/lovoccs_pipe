
# == This file is for storing filepaths to different folders (databases) which are used by lv_pipe == #
# WARNING: This is read by Bash; so no spaces when defining variables!!!


# == bash configurables for photometric_correction.sh == #

# a list of possible catalogs and the corresponding instruments
# this is only used during photometric_correction.sh... tweaking the refcats in LSSTPipe requires tweaking the loops and options in run_steps_Gen3.sh 
POSSIBLE_CATALOGS=("des_dr2" "legacy_survey_dr9" "ps1_dr1" "sm_dr1" "sm_dr2")
INSTRUMENTS=("des" "legacy" "ps1" "sm" "sm")

# filepath to refcats, assumes refcat for a given cluster is stored in ${CAT_DB}/${CLN}/...
CAT_DB='/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new'

# filepath to extinction; assumes fits w. extinction is stored in ${EXT_DB}/${CLN}.fits
EXT_DB='/gpfs/data/idellant/Clusters/galactic_extinction_database/data_extinction_irsa'


# == bash configurables for photo_z.sh == #

# filepath to specz; looks for csv's stored in ${SPECZ_DB}/${CLN}_ned_select.csv
SPECZ_DB='/gpfs/data/idellant/Clusters/spec_z_database'


# == bash configurables for mass_map.sh == #

#TODO change the structure of how we save xray data to make it friendlier for collaborators to use pipe
# filepath to xray db; currently searches for ${XRAY_DB}/${CLN}/Chandra/broad_flux_smoothed.fits
XRAY_DB='/gpfs/data/idellant/Clusters/Xray'





