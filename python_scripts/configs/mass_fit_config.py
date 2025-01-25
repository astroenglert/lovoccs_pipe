
import numpy as np

# == config for mass_fit.py == #

# resolution of different instruments,
instrument_resolution = {
                         'decam' : 0.263,
                         'hsc' : 0.168,
                        }

# quality cuts to be applied...
#TODO re-write this and the corresponding function in mass_fit to select from an interval rather than parse an inequality
quality_cuts = [
                ('r_cmodel_magerr','<',(np.log(10)/2.5)/5), # SN-cut written a little weird since dm ~ df/f
                ('blendedness','<',0.42),
                ('res','>',0.3),
                ('sigmae','<',0.4),
                ('e1','<',4),
                ('e2','<',4), # these two together are equivalent to |e| < 4
                ('g1','<',2),
                ('g2','<',2), # and these two together are equivalent to |g| < 2
                ('mod_chi2','<',4), # beware, in our old catalogs this is named chi2_mod, new is mod_chi2
                ('odds','>',0.95),
                ('z_phot','>',0.15), #TODO should this be cluster-dependant?
                ('z_phot','<',1.4),
                ('r_cmodel_mag','<',26),
                ('r_cmodel_mag','>',17),
               ]


# options for bootstrapping the peak px in mass_map
# bootstrap the peak position?
bootstrap_peak = True

# filter options for mass_map calls
map_filter = 'schirmer'



