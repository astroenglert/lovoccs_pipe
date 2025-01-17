# TIO -> Tools for IO (we love abreviations in astronomy)

# return the headers for a specific catalog
# eventually we should rewrite this into something which loads info from a json or smthng?
def get_instrument_headers(catalog_name):
    '''
    This is a temporary function which outputs a hand-coded dictionary of column names, eventually we need to replace this with a full solution which loads the appropriate configs from jsons or some other format.
    
    Args:
      catalog_name: String; one of decam, des, sm, ps1, sdss from which to load the proper headers for a reference catalog
    
    Returns:
        output_dict: Dictionary; a dictionary which maps ra/dec and magnitudes to the appropriate column of a reference catalog
    
    '''

    #output_dict = None
    if catalog_name == "decam":
        output_dict = {'ra_name':'ra',
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
                       }
    
    elif catalog_name == "des":
        output_dict = {'ra_name':'ra',
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
                       }
    
    elif catalog_name == "sm":
        output_dict = {'ra_name':'raj2000',
                       'dec_name':'dej2000',
                       'u_psf_mag':'v_psf',
                       'g_psf_mag':'g_psf',
                       'r_psf_mag':'r_psf',
                       'i_psf_mag':'i_psf',
                       'z_psf_mag':'z_psf',
                       'Y_psf_mag':None,
                       'u_cmd_mag':None,
                       'g_cmd_mag':None,
                       'r_cmd_mag':None,
                       'i_cmd_mag':None,
                       'z_cmd_mag':None,
                       'Y_cmd_mag':None,
                       }
    
    elif catalog_name == "ps1":
        output_dict = {'ra_name':'RAJ2000',
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
                       }
    
    elif catalog_name == "sdss":
        output_dict = {'ra_name':'RAJ2000',
                       'dec_name':'DEJ2000',
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
                       }
    
    else:
        raise Exception(f'Dictionary not found for {catalog_name}')
    
    return output_dict

