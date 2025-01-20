
# == This file is for storing filepaths to different folders (databases) which are used by lv_pipe == #
# WARNING: This is read by Bash; so no spaces when defining variables!!!

# == filepaths for photometric_correction.sh == #

# filepath to refcats, assumes refcat for a given cluster is stored in ${CAT_DB}/${CLN}/...
CAT_DB='/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new'

# filepath to extinction; assumes fits w. extinction is stored in ${EXT_DB}/${CLN}.fits
EXT_DB='/gpfs/data/idellant/Clusters/galactic_extinction_database/data_extinction_irsa'

