# calibrate w. metadetect

import sys
import os
from pathlib import Path

from multiprocessing import Pool

import numpy as np

from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.stats import binned_statistic

import matplotlib.pyplot as pl

from astropy.io import ascii
from astropy.table import Table, hstack, vstack

# Homebrew modules here
from ..mass_map.mass_map import load_quality_cuts

# loading config
from ..configs.meta_calibrate_config import quality_cuts, shape_type, delta_gamma


#TODO define a version of this function to allow for weights
# let's define a function to compute the responsivity 
def get_responsivity_unweighted(tab_1p,tab_1m,tab_2p,tab_2m,shape_type='sdss',delta_gamma=0.02,bin_resp=3,resp_min=0.1):
    '''
    A helper function to compute the responsivity given metadetect catalogs
    
    Args:
        tab_1p: Astropy Table; table from positive-e1 catalogs
        tab_1m: Astropy Table; table from negative-e1 catalogs
        tab_2p: Astropy Table; table from positive-e2 catalogs
        tab_2m: Astropy Table; table from negative-e2 catalogs
        shape_type: String; type of shape measurement to use for estimating g1/g2
        delta_gamma: Float; Change in shear between positive/negative catalogs
        bin_resp: Int; number of bins in resolvedness used to compute resp
        resp_min: Float; minimum responisivty to bin
        
    Returns:
        [R11,R12,R21,R22]: array of numpy arrays; contains the binned responsivities
        [R11_err,R12_err,R21_err,R22_err]: array of numpy arrays; contains the error on binned responsivities
    
    '''
    
    # use binned_statistic to take care of everything binning-related :)
    # bins have to be identical here
    bins = np.linspace(resp_min,1,bin_resp+1)
    e1_1p,_,_ = binned_statistic(values=tab_1p[shape_type + '_e1'],x=tab_1p[shape_type + '_res'],bins=bins,statistic='mean')
    e2_1p,_,_ = binned_statistic(values=tab_1p[shape_type + '_e2'],x=tab_1p[shape_type + '_res'],bins=bins)
    e1_1m,_,_ = binned_statistic(values=tab_1m[shape_type + '_e1'],x=tab_1m[shape_type + '_res'],bins=bins)
    e2_1m,_,_ = binned_statistic(values=tab_1m[shape_type + '_e2'],x=tab_1m[shape_type + '_res'],bins=bins)
    e1_2p,_,_ = binned_statistic(values=tab_2p[shape_type + '_e1'],x=tab_2p[shape_type + '_res'],bins=bins)
    e2_2p,_,_ = binned_statistic(values=tab_2p[shape_type + '_e2'],x=tab_2p[shape_type + '_res'],bins=bins)
    e1_2m,_,_ = binned_statistic(values=tab_2m[shape_type + '_e1'],x=tab_2m[shape_type + '_res'],bins=bins)
    e2_2m,_,_ = binned_statistic(values=tab_2m[shape_type + '_e2'],x=tab_2m[shape_type + '_res'],bins=bins)
    
    # helper function to get a standard error
    stderr = lambda x : np.std(x)/np.sqrt(len(x))
    
    # also measure the error on these means to estimate the error on each resp element
    e1_1p_err,_,_ = binned_statistic(values=tab_1p[shape_type + '_e1'],x=tab_1p[shape_type + '_res'],bins=bins,statistic=stderr)
    e2_1p_err,_,_ = binned_statistic(values=tab_1p[shape_type + '_e2'],x=tab_1p[shape_type + '_res'],bins=bins,statistic=stderr)
    e1_1m_err,_,_ = binned_statistic(values=tab_1m[shape_type + '_e1'],x=tab_1m[shape_type + '_res'],bins=bins,statistic=stderr)
    e2_1m_err,_,_ = binned_statistic(values=tab_1m[shape_type + '_e2'],x=tab_1m[shape_type + '_res'],bins=bins,statistic=stderr)
    e1_2p_err,_,_ = binned_statistic(values=tab_2p[shape_type + '_e1'],x=tab_2p[shape_type + '_res'],bins=bins,statistic=stderr)
    e2_2p_err,_,_ = binned_statistic(values=tab_2p[shape_type + '_e2'],x=tab_2p[shape_type + '_res'],bins=bins,statistic=stderr)
    e1_2m_err,_,_ = binned_statistic(values=tab_2m[shape_type + '_e1'],x=tab_2m[shape_type + '_res'],bins=bins,statistic=stderr)
    e2_2m_err,_,_ = binned_statistic(values=tab_2m[shape_type + '_e2'],x=tab_2m[shape_type + '_res'],bins=bins,statistic=stderr)
    
    R11 = (e1_1p - e1_1m)/delta_gamma
    R12 = (e1_2p - e1_2m)/delta_gamma
    R21 = (e2_1p - e2_1m)/delta_gamma
    R22 = (e2_2p - e2_2m)/delta_gamma
    
    R11_err = np.sqrt( e1_1p_err**2 + e1_1m_err**2 )/delta_gamma
    R12_err = np.sqrt( e1_2p_err**2 + e1_2m_err**2 )/delta_gamma
    R21_err = np.sqrt( e2_1p_err**2 + e2_1m_err**2 )/delta_gamma
    R22_err = np.sqrt( e2_2p_err**2 + e2_2m_err**2 )/delta_gamma
    
    return [R11,R12,R21,R22], [R11_err,R12_err,R21_err,R22_err]


# and define an additional function to actually calibrate the ellipticities
def calibrate_shapes(tab_noshear,resp,diagonal=True,shape_type='sdss',bin_resp=3,min_res=0.1):
    '''
    Calibrate shapes using the responsivity
    
    Args:
        tab_noshear: Astropy Table; the noshear metadetect catalog to calibrate
        resp: numpy array; An array containing the responsivity [R11,R12,R21,R22]
        shape_type: String; the type of shape to use for calibration
        diagonal: Bool; Assume the responsivity is diagonal (usually true...)
        bin_resp: Int; number of bins in resolvedness used to compute resp
    
    Returns:
        g1: array-like; first comp of the calibrated reduced shear
        g2: array-like; second comp of the calibrated reduced shear
    
    '''
    
    # very simple algorithm here, just compute g1/g2 based on the resolvedness bins
    e1 = tab_noshear[shape_type + '_e1']
    e2 = tab_noshear[shape_type + '_e2']
    res = tab_noshear[shape_type + '_res']
    
    # binning this is a little annoying for an arbitrary number of bins, but not too horrible
    bins = np.linspace(min_res,1,bin_resp+1)
    select = []
    calib11 = []
    calib12 = []
    calib21 = []
    calib22 = []
    
    for i in range(bin_resp):
        
        # first setup the selection criteria
        if i == 0:
            select.append(res > bins[0])
        else:
            select.append( (res > bins[i]) & (res < bins[i+1]) )
        
        if diagonal:
            # element-wise division
            calib11.append(lambda x : x / resp[0][i])
            calib22.append(lambda x : x / resp[3][i])
        else:
            # multiply by resp, then divide by the determinant in the appropriate combinations
            det = resp[0][i]*resp[3][i] - resp[1][i]*resp[2][i]
            calib11.append(lambda x : x * resp[0][i]/det)
            calib12.append(lambda x : x * resp[1][i]/det)
            calib21.append(lambda x : x * resp[2][i]/det)
            calib22.append(lambda x : x * resp[3][i]/det)
    
    if diagonal:
        g1 = np.piecewise(e1,select,calib11)
        g2 = np.piecewise(e2,select,calib22)
    else:
        g1 = np.piecewise(e1,select,calib22) - np.piecewise(e2,select,calib12)
        g2 = - np.piecewise(e1,select,calib21) + np.piecewise(e2,select,calib11)
    
    return g1,g2


# draw quality check plots to verify that everything is okay
def quality_check(tab,output_directory,shape_type='sdss',min_res=0.1,max_blend=0.1):
    """
    
    Args:
      tab: Astropy Table; table to take shapes from
      output_directory: String; string pointing to the directory to write files
      shape_type: String; shape type used
      min_res: float; minimum resolvedness used
      max_blend: float; maximum blendedness used
    
    Returns:
      None
    
    """
    
    # load shears and other statistics
    g1 = tab['g1']
    g2 = tab['g2']
    res = tab[shape_type + '_res']
    sn = (2.5/np.log(10)) * (1/tab['r_cmodel_magerr'])
    blend = tab['blendedness']
    
    # define a helper function to get stderr
    stderr = lambda x : np.std(x)/np.sqrt(len(x))
    
    # bin and draw the quality plots!
    # first lets do resolvedness
    fig,ax = pl.subplots(2)
    
    res_bins = np.linspace(min_res,1,20)
    mean, edges, num = binned_statistic(x=res,values=g1,bins=res_bins)
    std, edges, num = binned_statistic(x=res,values=g1,bins=res_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[0].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    ax[0].set_ylim((-0.02,0.02))
    ax[0].set_ylabel("$ \overline{e_1} $")
    
    mean, edges, num = binned_statistic(x=res,values=g2,bins=res_bins)
    std, edges, num = binned_statistic(x=res,values=g2,bins=res_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[1].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    ax[1].set_ylim((-0.02,0.02))
    ax[1].set_ylabel("$ \overline{e_2} $")
    
    ax[0].set_title("Shear v. Resolvedness")
    ax[1].set_xlabel('Resolvedness')
    fig.savefig(output_directory + 'shear_v_res.png',bbox_inches='tight',dpi=720)
    
    # next let's do the SN
    fig,ax = pl.subplots(2)
    
    sn_bins = np.logspace(1,3,20)
    mean, edges, num = binned_statistic(x=sn,values=g1,bins=sn_bins)
    std, edges, num = binned_statistic(x=sn,values=g1,bins=sn_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[0].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    ax[0].set_ylim((-0.02,0.02))
    ax[0].set_ylabel("$ \overline{e_1} $")
    ax[0].set_xscale('log')
    
    mean, edges, num = binned_statistic(x=sn,values=g2,bins=sn_bins)
    std, edges, num = binned_statistic(x=sn,values=g2,bins=sn_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[1].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    ax[1].set_ylim((-0.02,0.02))
    ax[1].set_ylabel("$ \overline{e_2} $")
    ax[1].set_xscale('log')
    
    ax[0].set_title("Shear v. SN")
    ax[1].set_xlabel('SN')
    fig.savefig(output_directory + 'shear_v_sn.png',bbox_inches='tight',dpi=720)
    
    # and finally do this w. respect to the blendedness
    fig,ax = pl.subplots(2)
    
    blend_bins = np.linspace(0,max_blend,10)
    mean, edges, num = binned_statistic(x=blend,values=g1,bins=blend_bins)
    std, edges, num = binned_statistic(x=blend,values=g1,bins=blend_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[0].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    ax[0].set_ylim((-0.02,0.02))
    ax[0].set_ylabel("$ \overline{e_1} $")
    
    mean, edges, num = binned_statistic(x=blend,values=g2,bins=blend_bins)
    std, edges, num = binned_statistic(x=blend,values=g2,bins=blend_bins,statistic=stderr)
    edges = (edges[1:] + edges[:-1])/2
    
    ax[1].errorbar(x=edges,y=mean,yerr=std,capsize=3,fmt='.')
    ax[1].set_ylim((-0.02,0.02))
    ax[1].set_ylabel("$ \overline{e_2} $")
    
    ax[0].set_title("Shear v. Blendedness")
    ax[1].set_xlabel('Blendedness')
    fig.savefig(output_directory + 'shear_v_blend.png',bbox_inches='tight',dpi=720)
    
    return

if __name__ == '__main__':
    
    if len(sys.argv) == 7:
        
        file_noshear = sys.argv[1]
        file_1p = sys.argv[2]
        file_1m = sys.argv[3]
        file_2p = sys.argv[4]
        file_2m = sys.argv[5]
        output_directory = sys.argv[6]
        files = {
                 'noshear':file_noshear,
                 '1p':file_1p,
                 '1m':file_1m,
                 '2p':file_2p,
                 '2m':file_2m,
                }
    
    else:
        
        print("python meta_calibrate.py noshear 1p 1m 2p 2m output_directory")
        raise Exception("Improper Usage! Correct usage: python meta_calibrate.py noshear 1p 1m 2p 2m output_directory")
    
    # load in the catalogs and apply quality cuts 
    catalogs = {}
    for key,name in files.items():
        
        tab = ascii.read(name)
        tab = load_quality_cuts(tab,quality_cuts=quality_cuts)
        catalogs[key] = tab
    
    # compute the repsonsivity
    # Try by default to do 3 bins, but if the statistics don't support it decrease the bin number until Rii/Rii_err > 5 for all bins (or only one bin is left)
    bin_retry = True
    bin_resp = 3
    
    while bin_retry:
        
        # compute responsivity with the appropriate number of bins, then collect the "signal"
        resp, resp_err = get_responsivity_unweighted(catalogs['1p'],catalogs['1m'],catalogs['2p'],catalogs['2m'],shape_type=shape_type,delta_gamma=delta_gamma,bin_resp=bin_resp)
        signal11 = resp[0]/resp_err[0]
        signal22 = resp[3]/resp_err[3]
        
        # check that the "signal" in each responsivity bin is greater than 5
        if ( (np.sum(signal11 < 5) + np.sum(signal22 < 5)) == 0) or (bin_resp == 1):
            # only reaches here if it finished the loop! This is the number of bins it'll use
            bin_retry=False
            print(f'Binned into {bin_resp} bins.')
            print('R11, R12, R21, R22')
            print(resp)
            break
        else:
            # one of the "signals" was less than five, decrease the number of bins and retry!
            bin_resp = bin_resp - 1
            continue
    
    # compute the reduced shears!
    g1,g2 = calibrate_shapes(catalogs['noshear'],resp,shape_type=shape_type,bin_resp=bin_resp)
    
    tab = catalogs['noshear']
    
    # for the time being, do not apply any per-object weights
    tab['g1'] = g1
    tab['g2'] = g2
    tab['shape_weight'] = np.ones(len(tab))/(0.365**2)
    
    quality_check(tab,output_directory=output_directory)
    
    # and finally save this table!
    basename = Path(files['noshear']).stem
    tab.write(output_directory + basename + '_meta_cal.csv',format='ascii.csv',overwrite=True)
    
    
    
    


