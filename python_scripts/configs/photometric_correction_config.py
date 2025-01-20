
# == configuration for extinction_correction.py == #

# dictionary of factors to convert E(B-V) to the per-band extinction for various instruments
# we use IRSA by default, although there are other potentially better options we intend to explore...
# lastly this is keyed by instrument and then by a magnitude a filter is mapped to (according to the instrument_header dictionary)
ebv_dictionary = {
                  'decam' : {'u_psf_mag':3.86,'g_psf_mag':3.11,'r_psf_mag':2.09,'i_psf_mag':1.53,'z_psf_mag':1.17,'Y_psf_mag':1.02,'u_cmd_mag':3.86,'g_cmd_mag':3.11,'r_cmd_mag':2.09,'i_cmd_mag':1.53,'z_cmd_mag':1.17,'Y_cmd_mag':1.02},
                  'des' : {'u_psf_mag':3.86,'g_psf_mag':3.11,'r_psf_mag':2.09,'i_psf_mag':1.53,'z_psf_mag':1.17,'Y_psf_mag':1.02,'u_cmd_mag':3.86,'g_cmd_mag':3.11,'r_cmd_mag':2.09,'i_cmd_mag':1.53,'z_cmd_mag':1.17,'Y_cmd_mag':1.02},
                  'legacy' : {'u_psf_mag':3.86,'g_psf_mag':3.11,'r_psf_mag':2.09,'i_psf_mag':1.53,'z_psf_mag':1.17,'Y_psf_mag':1.02,'u_cmd_mag':3.86,'g_cmd_mag':3.11,'r_cmd_mag':2.09,'i_cmd_mag':1.53,'z_cmd_mag':1.17,'Y_cmd_mag':1.02},
                  'ps1' : {'u_psf_mag':None,'g_psf_mag':3.08,'r_psf_mag':2.20,'i_psf_mag':1.63,'z_psf_mag':1.28,'Y_psf_mag':1.07},
                  'sm' : {'u_psf_mag':3.92,'g_psf_mag':2.89,'r_psf_mag':2.22,'i_psf_mag':1.55,'z_psf_mag':1.17,'Y_psf_mag':None},
                  'sdss' : {'u_psf_mag':4.16,'g_psf_mag':3.20,'r_psf_mag':2.21,'i_psf_mag':1.64,'z_psf_mag':1.22,'Y_psf_mag':None},
                 }


# == configuration for color_terms.py == #

# dictionary mapping instrument headers to LSST magnitudes/columns
# I've used a shorthand here (cmd -> cmodel) just to keep the keys a little shorter
# adding additional refcats here just requires adding an additional entry to this dictionary
instrument_headers = {'decam' : {'ra_name':'ra',
                                 'dec_name':'dec',
                                 'u_psf_mag':'u_psf_mag',
                                 'g_psf_mag':'g_psf_mag',
                                 'r_psf_mag':'r_psf_mag',
                                 'i_psf_mag':'i_psf_mag',
                                 'z_psf_mag':'z_psf_mag',
                                 'Y_psf_mag':'Y_psf_mag',
                                 'u_cmd_mag':'u_cmodel_mag',
                                 'g_cmd_mag':'g_cmodel_mag',
                                 'r_cmd_mag':'r_cmodel_mag',
                                 'i_cmd_mag':'i_cmodel_mag',
                                 'z_cmd_mag':'z_cmodel_mag',
                                 'Y_cmd_mag':'Y_cmodel_mag',
                                },
                      'des' : {'ra_name':'ra',
                               'dec_name':'dec',
                               'u_psf_mag':None,
                               'g_psf_mag':'wavg_mag_psf_g',
                               'r_psf_mag':'wavg_mag_psf_r',
                               'i_psf_mag':'wavg_mag_psf_i',
                               'z_psf_mag':'wavg_mag_psf_z',
                               'Y_psf_mag':'wavg_mag_psf_y',
                               'u_cmd_mag':None,
                               'g_cmd_mag':None,
                               'r_cmd_mag':None,
                               'i_cmd_mag':None,
                               'z_cmd_mag':None,
                               'Y_cmd_mag':None,
                              },
                      'legacy' : {'ra_name':'ra',
                                  'dec_name':'dec',
                                  'u_psf_mag':None,
                                  'g_psf_mag':'mag_g',
                                  'r_psf_mag':'mag_r',
                                  'i_psf_mag':None,
                                  'z_psf_mag':'mag_z',
                                  'Y_psf_mag':None,
                                  'u_cmd_mag':None,
                                  'g_cmd_mag':None,
                                  'r_cmd_mag':None,
                                  'i_cmd_mag':None,
                                  'z_cmd_mag':None,
                                  'Y_cmd_mag':None,
                                 },
                      'sm' : {'ra_name':'raj2000',
                              'dec_name':'dej2000',
                              'u_psf_mag':'v_psf',
                              'g_psf_mag':'g_psf',
                              'r_psf_mag':'r_psf',
                              'i_psf_mag':'i_psf',
                              'z_psf_mag':'z_psf',
                              'u_cmd_mag':None,
                              'g_cmd_mag':None,
                              'r_cmd_mag':None,
                              'i_cmd_mag':None,
                              'z_cmd_mag':None,
                              'Y_cmd_mag':None,                                
                             },
                      'ps1' : {'ra_name':'RAJ2000',
                               'dec_name':'DEJ2000',
                               'u_psf_mag':None,
                               'g_psf_mag':'gmag',
                               'r_psf_mag':'rmag',
                               'i_psf_mag':'imag',
                               'z_psf_mag':'zmag',
                               'Y_psf_mag':'ymag',
                               'u_cmd_mag':None,
                               'g_cmd_mag':None,
                               'r_cmd_mag':None,
                               'i_cmd_mag':None,
                               'z_cmd_mag':None,
                               'Y_cmd_mag':None,
                              },
                      'sdss' : {'ra_name':'RA_ICRS',
                                'dec_name':'DE_ICRS',
                                'u_psf_mag':'upmag',
                                'g_psf_mag':'gpmag',
                                'r_psf_mag':'rpmag',
                                'i_psf_mag':'ipmag',
                                'z_psf_mag':'zpmag',
                                'Y_psf_mag':None,
                                'u_cmd_mag':None,
                                'g_cmd_mag':None,
                                'r_cmd_mag':None,
                                'i_cmd_mag':None,
                                'z_cmd_mag':None,
                                'Y_cmd_mag':None,
                               },
                     }

# similar to instrument headers, but this is specifically for the reference stars
std_star_headers = {'decam' : {'u_psf_mag':'u_DECam',
                                 'g_psf_mag':'g_DECam',
                                 'r_psf_mag':'r_DECam',
                                 'i_psf_mag':'i_DECam',
                                 'z_psf_mag':'z_DECam',
                                 'Y_psf_mag':'Y_DECam',
                                 'u_cmd_mag':None,
                                 'g_cmd_mag':None,
                                 'r_cmd_mag':None,
                                 'i_cmd_mag':None,
                                 'z_cmd_mag':None,
                                 'Y_cmd_mag':None,
                                },
                      'sm' : {
                              'u_psf_mag':'v_SM',
                              'g_psf_mag':'g_SM',
                              'r_psf_mag':'r_SM',
                              'i_psf_mag':'i_SM',
                              'z_psf_mag':'z_SM',
                              'u_cmd_mag':None,
                              'g_cmd_mag':None,
                              'r_cmd_mag':None,
                              'i_cmd_mag':None,
                              'z_cmd_mag':None,
                              'Y_cmd_mag':None,                                
                             },
                      'ps1' : {
                               'u_psf_mag':None,
                               'g_psf_mag':'g_PS1',
                               'r_psf_mag':'r_PS1',
                               'i_psf_mag':'i_PS1',
                               'z_psf_mag':'z_PS1',
                               'Y_psf_mag':'y_PS1',
                               'u_cmd_mag':None,
                               'g_cmd_mag':None,
                               'r_cmd_mag':None,
                               'i_cmd_mag':None,
                               'z_cmd_mag':None,
                               'Y_cmd_mag':None,
                              },
                      'sdss' : {
                                'u_psf_mag':'u_SDSS',
                                'g_psf_mag':'g_SDSS',
                                'r_psf_mag':'r_SDSS',
                                'i_psf_mag':'i_SDSS',
                                'z_psf_mag':'z_SDSS',
                                'Y_psf_mag':None,
                                'u_cmd_mag':None,
                                'g_cmd_mag':None,
                                'r_cmd_mag':None,
                                'i_cmd_mag':None,
                                'z_cmd_mag':None,
                                'Y_cmd_mag':None,
                               },
                     }

# a dictionary directing the code to colorterms stored on-disk
colorterm_dictionary = {
                        'des':'/oscar/data/idellant/Clusters/color_term/DES_factor.csv',
                        'legacy':'/oscar/data/idellant/Clusters/color_term/DES_factor.csv',
                        'ps1':'/oscar/data/idellant/Clusters/color_term/PS1_factor.csv',
                        'sm':'/oscar/data/idellant/Clusters/color_term/SM_factor.csv',
                        'sdss':'/oscar/data/idellant/Clusters/color_term/SDSS_factor.csv',
                       }


#TODO we should make things like the color-cuts used configurable by the user

# == configuration for zero_point.py == #

use_locus = ['u']







