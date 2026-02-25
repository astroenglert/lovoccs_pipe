
import numpy as np

# == configuration for mass_map.py == #

# resolution of different instruments,
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
                ('g1','<',2),
                ('g2','<',2), # and these two together are equivalent to |g| < 2
                ('mod_chi2','<',4), # beware, in our old catalogs this is named chi2_mod, new is mod_chi2
                ('odds','>',0.7),
                ('z_phot','>',0.15), #TODO should this be cluster-dependant?
                ('z_phot','<',1.4),
                ('r_cmodel_mag','<',26),
                ('r_cmodel_mag','>',17),
               ]

# range of aperture sizes to consider
aperture_sizes = np.arange(3000,30000,1000)

# keys are name for SourceCatalog object, values are what we'll name the columns
# sn_columns are measurements made (via photutils SourceCatalog) on the SN-maps
sn_columns = {
              'label':'SourceID',
              'max_value':'SN_peak',
              'maxval_xindex':'x_sn_max',
              'maxval_yindex':'y_sn_max',
              'xcentroid':'x_sn_centroid',
              'ycentroid':'y_sn_centroid',
              'covar_sigx2':'xx_sn_gauss',
              'covar_sigxy':'xy_sn_gauss',
              'covar_sigy2':'yy_sn_gauss',
              'area':'area',
              'orientation':'sn_theta',
              'fwhm':'sn_fwhm',
              'ellipticity':'sn_ellip',
             }

# map_columns are measurements made (via photutils SourceCatalog) on the SN-maps
map_columns = {
               'max_value':'Map_max',
               'maxval_xindex':'x_ap_max',
               'maxval_yindex':'y_ap_max',
               'xcentroid':'x_ap_centroid',
               'ycentroid':'y_ap_centroid',
               'covar_sigx2':'xx_ap_gauss',
               'covar_sigxy':'xy_ap_gauss',
               'covar_sigy2':'yy_ap_gauss',
               'orientation':'ap_theta',
               'fwhm':'ap_fwhm',
               'ellipticity':'ap_ellip',
              }


