
import numpy as np

# == config for red_sequence.py == #

# cut for SN
sn_cut = 10

# tolerance for selecting cluster-members
zspec_tolerance = 0.01

# filter map for specifying columns, and tag for errors
error_tag = 'err'
filter_map = {
              'u':'u_cmodel_mag',
              'g':'g_cmodel_mag',
              'r':'r_cmodel_mag',
              'i':'i_cmodel_mag',
              'z':'z_cmodel_mag',
             }


# instrument resolution dictionary
instrument_resolution = {
                         'decam' : 0.263,
                         'hsc' : 0.168,
                        }

# the color bins for gr/ri selection
color_bins_gr = np.linspace(0.65,1,7)
color_bins_ri = np.linspace(0.25,0.45,6)

# magnitude bins (script assumes these are linspaced)
mag_bins = np.linspace(13.5,19.5,21)

# scale magnitude based on the cluster redshift
scale_upper_mag = True
z_pedestal = 0.08

