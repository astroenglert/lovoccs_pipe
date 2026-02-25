
import sys

import numpy as np

import astropy.io.fits as fits
import healpy as hp
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord

# homebrew modules here
from .color_terms import get_instrument_headers
from ..configs.photometric_correction_config import ebv_dictionary

def get_extinction(ext_img,coords,catalog,ext_map):
    '''
    This loads an extinction from some database. For the time being, it's setup to work for our IRSA_EBV extinction maps, but we need to configure this further in the future
    
    Args:
        ext_img: String; filepath to the extinction image
        coords: Skycoord; list of coordinates from a catalog to apply the correction too
        catalog: Astropy Table; the catalog to correct
        ext_map: String; string specifying the correction to apply, currently not in use
    
    Returns:
        cut_catalog: A cut version of the catalog within the FOV of the extinction map ext_img
        ebv: The BV of every source in cut_catalog
        ebv_error: [opt] fetched error of BV
    
    '''
    
    # Written for IRSA extinction by default
    # To add more extinction maps, expand this function!
    if ext_map=='IRSA_EBV':
        ext=0
        with fits.open(ext_img) as hdul:
            wcs = WCS(hdul[ext].header)
            index = wcs.world_to_array_index(coords)
            
            # the FOV of the refcat could be larger than the extinction map
            # so we truncate the output catalog to the size of the IRSA image
            ncol = hdul[ext].header["NAXIS1"]
            nrow = hdul[ext].header["NAXIS2"]
            
            index0 = np.array( index[0] )
            index1 = np.array( index[1] )
            
            select =  index0 >= 0
            select &= index1 >= 0
            select &= index0 < nrow
            select &= index1 < ncol
            
            index = ( index0[select], index1[select] )
            
            ebv = hdul[ext].data[index]
            ebv_err = None

        cut_catalog = catalog[select]
    
    if ext_map=='PLANCK2016':
        ext=1
        with fits.open(ext_img) as hdul:
            data = hdul[ext].data  
            tau353 = data["TAU353"]
            err_tau = data["ERR_TAU"]
            ebv_map = tau353 * 1.49e4
            ebv_err_map = err_tau * 1.49e4
            nside = hp.npix2nside(len(ebv_map))
            
            theta = 0.5 * np.pi - coords.dec.radian   # colatitude
            phi = coords.ra.radian                    # longitude
            pix = hp.ang2pix(nside, theta, phi)

            # Interpolate E(B–V) at each coordinate
            ebv = hp.get_interp_val(ebv_map, theta, phi)
            # ebv_err = hp.get_interp_val(ebv_err_map, theta, phi)

            # propagating interpolation error with ebv_err
            ebv_var_map = ebv_err_map**2
            tiles, weights = hp.get_interp_weights(nside,theta,phi=phi)
            ebv_err = np.sqrt(np.sum(ebv_var_map[tiles] * weights**2, axis=0))
        cut_catalog = catalog.copy()

    # returns the cut catalog and EBV to be applied to each coordinate
    return cut_catalog,ebv,ebv_err


def load_extinction_transformation(instrument,band,ext_map):
    '''
    To apply an extinction correction, we have to transform the EBV into the correct term for a given band. Eventually we'll set this function up so that it loads from a proper config file, rather than the harcoded coefficients below.
    
    Args:
        instrument: String; specifis the instrument/refcat you're applying the correction too
        band: String; for now this is a string which specifies the band the correction is applied in
        ext_map: String; OPTIONAL; this specifies the ext_map. Not in use outside the default for now
    
    Returns:
        extinction_coefficient: Float; the correction to go from EBV extinction to extinction in the appropriate band
    
    '''
    
    if ext_map in ['IRSA_EBV','PLANCK2016']:
        extinction_coefficient = ebv_dictionary[instrument][band]
    else:
        raise Exception(f'No map {exp_map}')
    
    return extinction_coefficient
    
# apply the extinction correction to a catalog and return the dereddened version
def apply_extinction(catalog,catalog_name,ext_img,ext_map):
    '''
    Applies the extinction correction to catalog, 
    
    Args:
        catalog: Astropy Table; a catalog with headers defined in get_instrument_headers()
        catalog_name: String; the string specifying which of the instrument headers to use
        ext_map: String; specified the extinction to correct for, defaults to ISRA_EBV
    
    Returns:
        dereddened_catalog: Astropy Table; a copy of the original table which has been dereddened
    
    '''
    
    # make a copy of the catalog
    catalog_copy = catalog.copy()
    
    # load headers from the target catalog
    headers = get_instrument_headers(catalog_name)

    # get skycorr object containing it
    coords = SkyCoord(catalog_copy[headers['ra_name']],catalog_copy[headers['dec_name']],frame='icrs',unit='deg')
    
    # get the extinction correction for this cluster
    cut_catalog, ebv, ebv_err = get_extinction(ext_img,coords,catalog_copy,ext_map)
    
    # create a list of the mags which exist for this catalog
    keys = list(headers.keys())
    mags = []
    for entry in keys:
        if entry[-4:] == '_mag':
            mags.append(entry)
    
    # iterate through bands which aren't mapped to NONE
    # apply the extinction correction to each band with the correction from load_ext_trnsfm
    for mag in mags:
        
        # get the header for this catalog, if None continue
        magh = headers[mag]
        
        if magh == None:
            continue
        
        trans = load_extinction_transformation(catalog_name,mag,ext_map)
        
        # if this mag is in the input table, apply the correction
        if magh in cut_catalog.colnames:
            errh = magh + 'err'
            cut_catalog[magh] -= (ebv*trans)
            if errh in cut_catalog.colnames:
                cut_catalog[errh] = np.sqrt(cut_catalog[errh]**2 + ebv_err**2*trans**2)

    # return the dereddened catalog
    return cut_catalog
    
if __name__=='__main__':

    if len(sys.argv)!=5:
        print("python this.py ebv_image_filename catalog_csv_input_filename catalog_csv_output_filename instrument")
        raise Exception("Improper usage! Correct usage is python this.py ebv_image_filename catalog_csv_input_filename catalog_csv_output_filename instrument")
    
    # collecting arguments from cln
    ext_img = sys.argv[1]
    catalog = sys.argv[2]
    catalog_name = sys.argv[4]
    dereddened_filename = sys.argv[3]
    catalog = Table.read(catalog,format='ascii.csv')
    ext_map_choice = 'PLANCK2016'
    
    ext_img = '/gpfs/data/idellant/Clusters/dustmaps/planck/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
    # later changed in processing_step_configs.sh

    # deredding
    dered_cat = apply_extinction(catalog,catalog_name,ext_img,ext_map=ext_map_choice)
    dered_cat.write(dereddened_filename, format="ascii.csv", overwrite=True)

