
# == configuration for extinction_correction.py == #

# filter_map mapping template filters to columns from LSSTPipe
filter_map = {
              'u':'u_cmodel_mag',
              'g':'g_cmodel_mag',
              'r':'r_cmodel_mag',
              'i':'i_cmodel_mag',
              'z':'z_cmodel_mag',
             }

# tag applied to the end of LSSTPipe columns which specify errors
error_tag = 'err'

# location of cached fluxes to load; if None loads from module resources
#WARNING: the cached fluxes are DECam-fluxes, so if using a different instrument you should recompute them
external_cache = None

# used for cutting the matched catalog of specz-photoz when drawing plots and computing statistics
mod_chi2_cut = 4
odds_cut = 0.95

# 'ngvs' -> NGVS photo-z priors; 'cosmos' -> COSMOS photo-z priors
# more can be implemented by editing python_scripts/photo_z/sed_template.py
prior_choice = 'ngvs'


