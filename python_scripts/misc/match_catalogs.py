# Script to match catalogs together

# Eventually we should do something a little smarter...
# For now it's just matching to w/in 0.5" default, I'll also copy the temporary catalog-info function

import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as pl

from astropy.table import Table, hstack, unique
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

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
        output_dict = {'ra_name':'RA_ICRS',
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
                       }
    elif catalog_name == "legacy":
        output_dict = {'ra_name':'ra',
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
                       }
    else:
        raise Exception(f'Dictionary not found for {catalog_name}')
    
    return output_dict


def match_catalogs(refcat,refcat_inst,catalog,catalog_inst,refcat_tag='_ref',catalog_tag='_cat',sep=0.5*u.arcsec):
    '''
    Function which matches refcat to catalog; nearest neighbor within 0.5"
    
    Args:
        refcat: Astropy Table; reference catalog to match to
        refcat_inst: Dictionary; standard config dictionary specifying ra/dec headers
        catalog: Astronomy Table; catalog to be matched
        catalog_inst: Dictionary; standard config dictionary specifying ra/dec headers
        sep: Astropy Quantity; separation to match sources
    
    Returns:
        matched_catalog: Astropy Table; matched table
    '''
    
    # collect the headers
    refcat_headers = get_instrument_headers(refcat_inst)
    catalog_headers = get_instrument_headers(catalog_inst)
    print(refcat_headers)
    print(refcat_inst)
    print(refcat)
    
    # use SkyCoords for efficient matching
    refcat_coords = SkyCoord(ra=refcat[refcat_headers['ra_name']],dec=refcat[refcat_headers['dec_name']],unit='deg',frame='icrs')
    catalog_coords = SkyCoord(ra=catalog[catalog_headers['ra_name']],dec=catalog[catalog_headers['dec_name']],unit='deg',frame='icrs')
    
    # the actual matching (standard astropy shenanigans)
    # using match_coordinates_sky over match_to_catalog_sky since the former by default matches the nearest-neighbor
    idx,d2d,d3d = match_coordinates_sky(refcat_coords,catalog_coords)
    matched = d2d < sep
    refcat_matched = refcat[matched]
    catalog_matched = catalog[idx[matched]]
    
    # append tags to the headers from the matched catalogs
    for i in refcat_matched.colnames:
        refcat_matched[i].name = i + refcat_tag
    for i in catalog_matched.colnames:
        catalog_matched[i].name = i + catalog_tag
    
    # assembling the combined catalog
    combined_catalog = hstack([ catalog_matched, refcat_matched])
    
    return combined_catalog
    
def matched_catalog_histogram(refcat,refcat_inst,catalog,catalog_inst,sep=0.5*u.arcsec,figure_tag=''):
    ''''
    Function to draw histogram of offsets in ra/dec between matched objects.
    
    This isn't a great way of implementing this since it recalculates the matching to make the plot
    But, for us, anything that takes less than an hour to do is moderately fast so it's fine!
    
    Args:
        refcat: Astropy Table; reference catalog to match to
        refcat_inst: Dictionary; standard config dictionary specifying ra/dec headers
        catalog: Astronomy Table; catalog to be matched
        catalog_inst: Dictionary; standard config dictionary specifying ra/dec headers
        sep: Astropy Quantity; separation to match sources
    
    Returns:
        True if figure saved successfully
    '''
    
    # collect the headers
    refcat_headers = get_instrument_headers(refcat_inst)
    catalog_headers = get_instrument_headers(catalog_inst)

    # use SkyCoords for efficient matching
    refcat_coords = SkyCoord(ra=refcat[refcat_headers['ra_name']],dec=refcat[refcat_headers['dec_name']],unit='deg',frame='icrs')
    catalog_coords = SkyCoord(ra=catalog[catalog_headers['ra_name']],dec=catalog[catalog_headers['dec_name']],unit='deg',frame='icrs')
    
    # the actual matching (standard astropy shenanigans)
    # using match_coordinates_sky over match_to_catalog_sky since the former by default matches the nearest-neighbor
    idx,d2d,d3d = match_coordinates_sky(refcat_coords,catalog_coords)
    matched = d2d < sep
    refcat_coords = refcat_coords[matched]
    catalog_coords = catalog_coords[idx[matched]]
    distances = d2d[matched].arcsec
    
    # collect differences in position
    coord_diff_ra = refcat_coords.ra - catalog_coords.ra
    coord_diff_ra = coord_diff_ra.arcsec
    coord_diff_dec = refcat_coords.dec - catalog_coords.dec
    coord_diff_dec = coord_diff_dec.arcsec
    
    pl.figure(1)
    pl.scatter(coord_diff_ra,coord_diff_dec,marker='.',alpha=0.4)
    pl.xlabel("$ \\Delta \\alpha $ (arcsec)")
    pl.ylabel("$ \\Delta \\delta $ (arcsec)")
    pl.title(f'Matching {catalog_inst} to {refcat_inst}') 
    pl.savefig(f'{catalog_inst}_matched_to_{refcat_inst}_residuals' + figure_tag + '.png',dpi=720)
    
    pl.figure(2)
    pl.hist(distances,bins='auto',histtype='step')
    pl.xlabel("$ \\Delta \\theta $ (arcsec)")
    pl.ylabel(" Frequency ")
    pl.title(f'Matching {catalog_inst} to {refcat_inst}')
    pl.savefig(f'{catalog_inst}_matched_to_{refcat_inst}_histogram' + figure_tag + '.png',dpi=720)
    
    return True

if __name__ == '__main__':
    
    # check that two catalogs and two instrument names are passed (and optionally a distance)
    # then run match_catalogs() to match them and save the result
    # producing the histogram is mostly as a utility for debugging
    
    if len(sys.argv)!=9:
        print("Usage: python this.py refcat_filename refcat_inst refcat_tag catalog_filename catalog_inst catalog_tag combine_filename max_sep_arcsec")
        sys.exit(1)
    
    # loading cln arguments
    refcat_file = sys.argv[1]
    refcat_inst = sys.argv[2]
    refcat_tag = sys.argv[3]
    catalog_file = sys.argv[4]
    catalog_inst = sys.argv[5]
    catalog_tag = sys.argv[6]
    combine_file = sys.argv[7]
    max_sep = float(sys.argv[8])*u.arcsec
    
    # loading tables
    refcat = Table.read(refcat_file,format='ascii.csv')
    catalog = Table.read(catalog_file,format='ascii.csv')
    
    # matching catalogs
    combined = match_catalogs(refcat,refcat_inst,catalog,catalog_inst,refcat_tag,catalog_tag,max_sep)
    
    # writing the merged catalog
    combined.write(combine_file,format='ascii.csv',overwrite=True)
    
    
