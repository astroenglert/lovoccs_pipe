import sys
import os

from multiprocessing import Pool

import numpy as np

from scipy.integrate import cumulative_trapezoid as cumtrapz

import matplotlib.pyplot as pl

from astropy.io import ascii
from astropy.table import Table, hstack, vstack

# Homebrew modules here
from .hsc.utilities import create_calibs, new_columns
from .hsc.gen_hsc_calibrations import fix_nan

# loading config
#from ..configs.shear_calibration_config import 

#TODO quality-check plots for verifying shears
def draw_quality_check():
    
    
    
    pass


# this wraps the hsc correction into a single function call, in-case we want to test it further\
#TODO not sure why this isn't working, something way down in griddata is breaking, someone should fix this
def apply_hsc_correction(table):
    '''
    This is apples the hsc sbear correction algorithm (built for HSC-SSP). In our testing this doesn't work great for DECam, but it got us through the first two papers and worked well-enough!
    
    Args:
        table: Astropy Table; a table containing columns exported from LSSTPipe, including r_cmodel_magerr, res, and ei_psf_sdss
    
    Returns:
        output_table: Astropy Table; a table containing information from the shear calibration and the reduced-shears themselves.
    '''

    #TODO overhaul the hard-coded column-names for standard-IO
    shear_information = create_calibs(table)

    # Fix NaN/inf values, which could result from some NaN/inf values in the input catalog.
    # Typically there are very few of these.
    for entry in new_columns:
        fix_nan(shear_information, entry)
    
    e1 = shear_information["e1"] 
    e2 = shear_information["e2"]
    e_rms = shear_information["rms_e_d"]
    
    m = shear_information["shear_m"]
    c1 = shear_information["shear_c1"]
    c2 = shear_information["shear_c2"]
    
    weight = shear_information["shape_weight"]
    
    R = 1. - np.sum(weight * e_rms ** 2.)/np.sum(weight)
    
    m_mean = np.sum(weight * m)/np.sum(weight)
    c1_mean = np.sum(weight * c1)/np.sum(weight)
    c2_mean = np.sum(weight * c2)/np.sum(weight)
    print("R, m_mean, c1_mean, c2_mean: ", R, m_mean, c1_mean, c2_mean)
    
    g1_0 = ( e1 / (2. * R) - c1_mean ) / (1. + m_mean)
    g2_0 = ( e2 / (2. * R) - c2_mean ) / (1. + m_mean)
    
    # Another method
    # The difference is small because ci is small
    g1 = ( e1 / (2. * R) - c1 ) / (1. + m_mean)
    g2 = ( e2 / (2. * R) - c2 ) / (1. + m_mean)

    shear_information['g1'] = g1
    shear_information['g2'] = g2
    
    shear_information['g1_0'] = g1_0
    shear_information['g2_0'] = g2_0
    
    return shear_information


# this carries out a simple correction assuming a fixed responsivity (dell'Antonio+2020, Utsumi+2020)
def naive_responsivity_correction(table,c=0,m=0):
    '''
    This applies a naive 'shear_calibration' with options to include a fixed additive and multiplicative-bias. Assumes a fixed responsivity per-source consistent with dell'Antonio+2020 and Utsumi+2020.
    
    Args:
        table: Astropy Table; a table containing columns exported from LSSTPipe, including r_cmodel_magerr, res, and ei_psf_sdss
        c: float; an additive bias, defaults to zero (no-correction)
        m: float; a multiplicative bias, defaults to zero (no-correction)
    
    Returns:
        output_table: Astropy Table; a table containing the reduced shears g1/g2.
    '''
    output_table = Table()
    
    e1 = table['e1']
    e2 = table['e2']
    
    #TODO eventually we should compute this for our dataset properly, this is just a placeholder
    responsivity = 1-(0.365)**2
    
    g1 = (e1/(2*responsivity))
    g2 = (e2/(2*responsivity))
    
    output_table['g1'] = (g1 - c)/(1+m)
    output_table['g2'] = (g2 - c)/(1+m)
    
    return output_table

# a dictionary storing a list of the defined calibration methods
calibration_methods = {'hsc':apply_hsc_correction,'naive':naive_responsivity_correction}

if __name__ == '__main__':
    
    # check cln for arguments
    if len(sys.argv)!=4:
        print("python calibrate_shears.py catalog output_catalog_filename calibration_method")
        raise Exception("Improper Usage! Correct usage: python calibrate_shears.py catalog output_catalog_filename calibration_method")
        
    # collecting arguments from cln
    #TODO include an output-directory for figures?
    table_filename = sys.argv[1]
    output_catalog_filename = sys.argv[2]
    calib = sys.argv[3]
    
    if calib not in calibration_methods:
        print(f'WARNING: {calib} is not one of the allowed calibration methods, quitting the script now!')
        sys.exit()
    
    # load the table from disk
    table = ascii.read(table_filename)
    
    # use the appropriate shear-calibration algorithm to get an updated table
    shear_table = calibration_methods[calib](table)
    
    output_table = hstack([table,shear_table])
    
    # save the output table to disk
    output_table.write(output_catalog_filename,format='ascii.csv',overwrite=True)



    
    

