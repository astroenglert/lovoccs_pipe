
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
                ('blendedness','<',0.4),
                ('z_phot','>',0.1),
                ('z_phot','<',1.5),
                (f'{shape_type}_flag','==','False'), # catches any unphysical moments/shears
                (f'{shape_type}_res','>',0.3), # catches unresolved objects
                (f'{shape_type}_res','<',0.9), # weird class of objects have a large res but not good shears, cull these
               ]


# color-cuts to apply; this is exclusive to metadetect generally
# this applies cuts |color_cuts[0] - color_cuts[1]| < 5
color_cuts = [('g_cmodel_mag','r_cmodel_mag',10)]

