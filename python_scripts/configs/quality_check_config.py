
import numpy as np

# == config for quality_check.py == #

# mapping from ifilter to our filter-names (used for reading from query results)
filter_map = {
              'u': "u DECam c0006 3500.0 1000.0",
              'g': "g DECam SDSS c0001 4720.0 1520.0",
              'r': "r DECam SDSS c0002 6415.0 1480.0",
              'i': "i DECam SDSS c0003 7835.0 1470.0",
              'z': "z DECam SDSS c0004 9260.0 1520.0",
             }

# instrument resolution dictionary
instrument_resolution = {
                         'decam' : 0.263,
                         'hsc' : 0.168,
                        }

# quality cuts to be applied...
#TODO re-write this and the corresponding function in mass_map to select from an interval rather than parse an inequality
quality_cuts = [
                ('r_cmodel_magerr','<',(np.log(10)/2.5)/5), # SN-cut written a little weird since dm ~ df/f
                ('blendedness','<',0.42),
                ('res','>',0.3),
                ('sigmae','<',0.4),
                ('e1','<',4),
                ('e2','<',4), # these two together are equivalent to |e| < 4
               # ('g1','<',2),
               # ('g2','<',2), # and these two together are equivalent to |g| < 2
               # ('mod_chi2','<',4), # beware, in our old catalogs this is named chi2_mod, new is mod_chi2
               # ('odds','>',0.95),
               # ('z_phot','>',0.15), #TODO should this be cluster-dependant?
               # ('z_phot','<',1.4),
               # ('r_cmodel_mag','<',26),
               # ('r_cmodel_mag','>',17),
               ]

# list of bands to look for when generating quality-check plots
possible_bands = ['u','g','r','i','z']


