
# == config for meta_calibrate_config.py == #

import numpy as np

# metadetect has to use the same population for calculating R that is used for the eventual shears
# so selection-cuts are applied at this level rather than later in the process

#TODO make sure other shape algorithms are also compatible here
# shape measurement to use for metadetect
shape_type = 'sdss'

# change in shear between the +/- images
delta_gamma = 0.02

# Cuts to apply to the catalog
quality_cuts = [
                ('r_cmodel_magerr','<',(np.log(10)/2.5)/10), # SN-cut written a little weird since dm ~ df/f
                ('blendedness','<',0.1), # default 0.1
                ('mod_chi2','<',8), # weakened relative to defaults (now 4 -> 8) due to additional noise in metadetect
                ('odds','>',0.5), # weakened relative to defaults (now 0.7 -> 0.5) due to additional noise in metadetect
                ('z_phot','>',0.1),
                ('z_phot','<',1.5), # testing to remove A85's background cluster
                ('r_cmodel_mag','<',26),
                ('r_cmodel_mag','>',17),
                (f'{shape_type}_flag','==','False'), # catches any unphysical moments/shears
                (f'{shape_type}_res','>',0.1), # catches unresolved objects, defacto for meta is 0.1
               ]

