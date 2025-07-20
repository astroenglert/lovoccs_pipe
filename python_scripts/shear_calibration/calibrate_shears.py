import sys
import os

from multiprocessing import Pool

import numpy as np

from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.stats import binned_statistic

import matplotlib.pyplot as pl

from astropy.io import ascii
from astropy.table import Table, hstack, vstack

# Homebrew modules here
from .hsc.utilities import create_calibs, new_columns
from .hsc.gen_hsc_calibrations import fix_nan


# quality-check plots for verifying shears
def quality_check(tab,output_directory,min_res=0.3,max_blend=0.4):
    """
    
    Args:
      tab: Astropy Table; table to take shapes from
      output_directory: String; string pointing to the directory to write files
      min_res: float; minimum resolvedness used
      max_blend: float; maximum blendedness used
    
    Returns:
      None
    
    """
    
    # load shears and other statistics
    g1 = tab['g1']
    g2 = tab['g2']
    res = tab['res']
    sn = (2.5/np.log(10)) * (1/tab['r_cmodel_magerr'])
    blend = tab['blendedness']
    
    # define a helper function to get stderr
    stderr = lambda x : np.nanstd(x)/np.sqrt(np.sum(np.isfinite(x)))
    
    # bin and draw the quality plots!
    # first lets do resolvedness
    fig,ax = pl.subplots(2)
    
    res_bins = np.linspace(min_res,1,20)
    mean, edges, num = binned_statistic(x=res,values=g1,bins=res_bins,statistic=np.nanmean)
    std, edges, num = binned_statistic(x=res,values=g1,bins=res_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[0].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    #ax[0].set_ylim((-0.1,0.1))
    ax[0].set_ylabel("$ \overline{g_1} $")
    
    mean, edges, num = binned_statistic(x=res,values=g2,bins=res_bins,statistic=np.nanmean)
    std, edges, num = binned_statistic(x=res,values=g2,bins=res_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[1].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    #ax[1].set_ylim((-0.1,0.1))
    ax[1].set_ylabel("$ \overline{g_2} $")
    
    ax[0].set_title("Shear v. Resolvedness")
    ax[1].set_xlabel('Resolvedness')
    fig.savefig(output_directory + 'shear_v_res.png',bbox_inches='tight',dpi=720)
    
    # next let's do the SN
    fig,ax = pl.subplots(2)
    
    sn_bins = np.logspace(1,3,20)
    mean, edges, num = binned_statistic(x=sn,values=g1,bins=sn_bins,statistic=np.nanmean)
    std, edges, num = binned_statistic(x=sn,values=g1,bins=sn_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[0].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    #ax[0].set_ylim((-0.1,0.1))
    ax[0].set_ylabel("$ \overline{g_1} $")
    ax[0].set_xscale('log')
    
    mean, edges, num = binned_statistic(x=sn,values=g2,bins=sn_bins,statistic=np.nanmean)
    std, edges, num = binned_statistic(x=sn,values=g2,bins=sn_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[1].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    #ax[1].set_ylim((-0.1,0.1))
    ax[1].set_ylabel("$ \overline{g_2} $")
    ax[1].set_xscale('log')
    
    ax[0].set_title("Shear v. SN")
    ax[1].set_xlabel('SN')
    fig.savefig(output_directory + 'shear_v_sn.png',bbox_inches='tight',dpi=720)
    
    # and finally do this w. respect to the blendedness
    fig,ax = pl.subplots(2)
    
    blend_bins = np.linspace(0,max_blend,10)
    mean, edges, num = binned_statistic(x=blend,values=g1,bins=blend_bins,statistic=np.nanmean)
    std, edges, num = binned_statistic(x=blend,values=g1,bins=blend_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[0].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    #ax[0].set_ylim((-0.1,0.1))
    ax[0].set_ylabel("$ \overline{g_1} $")
    
    mean, edges, num = binned_statistic(x=blend,values=g2,bins=blend_bins,statistic=np.nanmean)
    std, edges, num = binned_statistic(x=blend,values=g2,bins=blend_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[1].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    #ax[1].set_ylim((-0.1,0.1))
    ax[1].set_ylabel("$ \overline{g_2} $")
    
    ax[0].set_title("Shear v. Blendedness")
    ax[1].set_xlabel('Blendedness')
    fig.savefig(output_directory + 'shear_v_blend.png',bbox_inches='tight',dpi=720)
    
    return



# this wraps the hsc correction into a single function call, in-case we want to test it further\
#TODO this keeps outputting nothing but zeros for resp/m/c; I might ask Soren to fix this since his build for SI/sledgehammer is working
def apply_hsc_correction(table):
    '''
    This is apples the hsc sbear correction algorithm (built for HSC-SSP). In our testing this doesn't work great for DECam, but it got us through the first two papers and worked well-enough!
    
    Args:
        table: Astropy Table; a table containing columns exported from LSSTPipe, including r_cmodel_magerr, res, and ei_psf_sdss
    
    Returns:
        output_table: Astropy Table; a table containing information from the shear calibration and the reduced-shears themselves.
    '''

    shear_information = create_calibs(table)

    # Fix NaN/inf values, which could result from some NaN/inf values in the input catalog.
    # Typically there are very few of these.
    for entry in new_columns:
        fix_nan(shear_information, entry)
    
    e1 = table["e1"] 
    e2 = table["e2"]
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
    print(g1)
    print(g2)
    
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
    
    # this is a naive calibration since it uses an established somewhat ad-hoc value for the resp
    # leave it to metadetect/hsc-calib to do this propery
    responsivity = 1-(0.365)**2
    
    g1 = (e1/(2*responsivity))
    g2 = (e2/(2*responsivity))
    
    output_table['g1'] = (g1 - c)/(1+m)
    output_table['g2'] = (g2 - c)/(1+m)
    output_table['shape_weight'] = np.ones(len(e1))/(0.365**2)
    
    return output_table

# a dictionary storing a list of the defined calibration methods
calibration_methods = {'hsc':apply_hsc_correction,'naive':naive_responsivity_correction}

if __name__ == '__main__':
    
    # check cln for arguments
    if len(sys.argv)!=5:
        print("python calibrate_shears.py catalog output_directory output_catalog_filename calibration_method")
        raise Exception("Improper Usage! Correct usage: python calibrate_shears.py catalog output_directory output_catalog_filename calibration_method")
        
    # collecting arguments from cln
    table_filename = sys.argv[1]
    output_directory = sys.argv[2]
    output_catalog_filename = sys.argv[3]
    calib = sys.argv[4]
    
    if calib not in calibration_methods:
        print(f'WARNING: {calib} is not one of the allowed calibration methods, quitting the script now!')
        sys.exit()
    
    # load the table from disk
    table = ascii.read(table_filename)
    
    # for hsc calibration, we need to run quality cuts a-priori to get a reliable calibration
    if calib == 'hsc':
        
        # unfortunately, I have to run cuts at this stage for HSC calibration to run properly
        # further cuts can be implemented in mass_map, but these are the bare-minimum needed to get a consistent calibration
        e = np.sqrt(table['e1']**2 + table['e2']**2)
        select = np.isfinite(e)
        select &= (e < 4)
        select &= (e > 0)
        select &= table['r_cmodel_mag'] < 26
        select &= table['r_cmodel_mag'] > 17
        select &= table["r_cmodel_magerr"] < 1.08574/5
        #select &= table['z_b'] > 0.15 # for checking a gen2 catalog!
        #select &= table['z_b'] > 0.15 # for checking a gen2 catalog!
        select &= table['z_phot'] > 0.15
        select &= table['z_phot'] < 1.4
        select &= table['res'] > 0.3
        select &= table['sigmae'] < 0.4
        select &= table['blendedness'] < 0.42
    
        table = table[select]
    
    # use the appropriate shear-calibration algorithm to get an updated table
    shear_table = calibration_methods[calib](table)
    
    output_table = hstack([table,shear_table])
    
    if calib == 'hsc':
        # run a final cut on g1/g2
        final_select = np.sqrt(output_table['g1']**2 + output_table['g2']**2) < 2
        output_table = output_table[final_select] 

    
    quality_check(output_table,output_directory)
    
    # save the output table to disk
    output_table.write(output_directory + output_catalog_filename,format='ascii.csv',overwrite=True)



    
    

