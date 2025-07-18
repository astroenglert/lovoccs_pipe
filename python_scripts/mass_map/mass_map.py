# run mass_map

import sys
import os
from pathlib import Path

from multiprocessing import Pool

import numpy as np

import matplotlib.pyplot as pl

import astropy.units as u
from astropy.io import ascii, fits
from astropy.table import Table, hstack, vstack
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.coordinates import SkyCoord

from astroquery.ipac.ned import Ned

from photutils.segmentation import detect_sources, SourceCatalog

from lsst.daf.butler import Butler

# HOMEBREW MODULES BELOW
from .filters import implemented_filters
from ..configs.mass_map_config import instrument_resolution, quality_cuts, aperture_sizes, sn_columns, map_columns
    

#TODO overhaul this to a standard-IO or function for running quality-cuts read from a config
def load_quality_cuts(table,quality_cuts=None):
    '''
    This is a temporary function which outputs a hand-coded dictionary of quality cuts we apply before passing data to mass_map. This requires no arguments, but in the future needs to be read from a config file
    
    Args:
        table: Astropy Table; a table to apply quality cuts to
        quality_cuts: array; an array of quality-cuts, stored in list. Formatting here is (COL,INEQU,CUT), e.g. (r_cmodel_mag,'<',10) is an instruction to preserve all sources with r_cmodel_mag below 10.
    
    Returns:
        cut_table: Astropy Table; a copy of table with the cuts applied

    '''
        
    # default to cuts specified in LVI
    if quality_cuts == None:
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
                        ('odds','>',0.95),
                        ('z_phot','>',0.15), #TODO should this be cluster-dependant?
                        ('z_phot','<',1.4),
                        ('r_cmodel_mag','<',26),
                        ('r_cmodel_mag','>',17),
                       ]
    
    select_sources = np.ones(len(table),dtype=bool)
    for cut in quality_cuts:
        # apply the appropriate inequality to select quality-sources
        if cut[1] == '<':  
            select_sources &= table[cut[0]] < cut[2]
        elif cut[1] == '>':
            select_sources &= table[cut[0]] > cut[2]
        elif cut[1] == '<=' or cut[1] == '=<':
            select_sources &= table[cut[0]] <= cut[2]
        elif cut[1] == '>=' or cut[1] == '=>':
            select_sources &= table[cut[0]] >= cut[2]
        elif cut[1] == '=' or cut[1] == '==':
            select_sources &= table[cut[0]] == cut[2]
        else:
            raise Exception('INEQ not recognized! Should be one of <,>,<=,>=,=')
    
    # and finally apply the cuts!  
    cut_table = table[select_sources]
    print("In total we have %s galaxies for mass-map"%(len(cut_table)))
    return cut_table
    

# compute the mass_map on a grid of points to be sampled
#TODO optimize this function, the nested for-loop is faster than numpy broadcasting (due to the sizes of the arrays it would have to create), but still should be a faster way of running this
#TODO right now we implicitly work in pixel-coordinates defined by the 12x12 skymap we define, this works well enough but will become an issue for studies beyond the central patches 

def compute_mass_map(x_grid,y_grid,x,y,g1,g2,weights,q_filter,filter_kwargs={}):
    '''
    This function computes the mass aperture-statistics at each point on a specified grid. Run quality-cuts, NaN filtering, etc. before this step!
    
    Args:
        x: Numpy array; an array of x-coordinates for each object
        y: Numpy array; an array of y-coordinates for each object
        x_grid: Numpy array; an NxM array of x-coordinates to sample the aperture-mass on
        y_grid: Numpy array; an NxM array of y-coordinates to sample the aperture-mass on
        g1; Numpy array; the shear g1 for each object
        g2; Numpy array; the shear g2 for each object
        weights: Numpy array; the weight for each object's shear
        q_filter; function; the filter-function used to compute Map
        kwargs; dict; kwargs passed to w_filter
    
    Returns:
        Map_E: Numpy array; an NxM array containing the E-mode aperture mass evaluated at each grid-point
        Map_B: Numpy array; an NxM array containing the B-mode aperture mass evaluated at each grid-point
        Map_V: Numpy array; an NxM array containing the variance in the aperture mass evaluated at each grid-point

    '''

    y_shape = len(y_grid[:,0])
    x_shape = len(x_grid[0,:])
    
    Map_E = np.zeros((y_shape,x_shape))
    Map_B = np.zeros((y_shape,x_shape))
    Map_V = np.zeros((y_shape,x_shape))
    
    if 'aperture_size' not in filter_kwargs:
        filter_area = np.pi * (8000)**2
    else:
        filter_area = np.pi * filter_kwargs['aperture_size']**2
    
    # an extra catch for an objects assigned NaN g1/g2
    #TODO should these be removed during shear calibration?
    nan_catch = np.isfinite(g1) & np.isfinite(g2)
    x = x[nan_catch]
    y = y[nan_catch]
    g1 = g1[nan_catch]
    g2 = g2[nan_catch]
    weights = weights[nan_catch]
    
    for i in range(y_shape):
        for j in range(x_shape):
            delta_x = x_grid[j,i] - x
            delta_y = y_grid[j,i] - y
            radius = np.sqrt(delta_x**2 + delta_y**2)
            theta = np.arctan2(delta_y,delta_x)
            g_T = -g1*np.cos(2*theta) - g2*np.sin(2*theta)
            g_X =  g1*np.sin(2*theta) - g2*np.cos(2*theta)
            g_mag = g1**2 + g2**2
            
            filter_values = q_filter(radius,**filter_kwargs)

            weight_sum = np.sum(weights)
            
            Map_E[i,j] = np.sum(filter_values*g_T*weights)*filter_area/weight_sum
            Map_B[i,j] = np.sum(filter_values*g_X*weights)*filter_area/weight_sum
            Map_V[i,j] = np.sum( (filter_values**2)*g_mag*(weights**2) )*(filter_area**2)/(2*(weight_sum**2))
    
    return Map_E, Map_B, Map_V


def compute_sn_px(x_sample,y_sample,x,y,g1,g2,weights,q_filter,filter_kwargs={}):
    '''
    Computes the SN-ratio at a specific pixel on the coadd
    
    Args:
        x_sample: float; the x-coordinate to evaluate Map at
        y_sample: float; the y-coordinate to evaluate Map at
        x; Numpy array; the shear x for each object
        y; Numpy array; the shear y for each object
        g1; Numpy array; the shear g1 for each object
        g2; Numpy array; the shear g2 for each object
        weights: Numpy array; the weight for each object's shear
        q_filter; function; the filter-function used to compute Map
        kwargs; dict; kwargs passed to w_filter
    
    Returns:
        sn: float; the signal-to-noise ratio at this pixel

    '''

    if 'aperture_size' not in filter_kwargs:
        filter_area = np.pi * (8000)**2
    else:
        filter_area = np.pi * filter_kwargs['aperture_size']**2
    
    # an extra catch for an objects assigned NaN g1/g2
    #TODO should these be removed during shear calibration?
    nan_catch = np.isfinite(g1) & np.isfinite(g2)
    x = x[nan_catch]
    y = y[nan_catch]
    g1 = g1[nan_catch]
    g2 = g2[nan_catch]
    weights = weights[nan_catch]
    

    delta_x = x_sample - x
    delta_y = y_sample - y
    radius = np.sqrt(delta_x**2 + delta_y**2)
    theta = np.arctan2(delta_y,delta_x)
    g_T = -g1*np.cos(2*theta) - g2*np.sin(2*theta)
    g_mag = g1**2 + g2**2
    
    filter_values = q_filter(radius,**filter_kwargs)
    
    weight_sum = np.sum(weights)
    
    Map_E = np.sum(filter_values*g_T*weights)*filter_area/weight_sum
    Map_V = np.sum( (filter_values**2)*g_mag*(weights**2) )*(filter_area**2)/(2*(weight_sum**2))
    sn = Map_E/np.sqrt(Map_V)
    
    return sn


#TODO For now, we're just using photutils but should we do something more sophisticated? Deblending for sub-peaks?
def detect_mass_peaks(Map_E,Map_B,Map_V,kwargs={'threshold':3,'npixels':25}):
    '''
    This function searches for peaks in the mass_map and creates a catalog of (arguably) useful statistics
    
    Args:
        Map_E: Numpy array; an NxM array containing the E-mode aperture mass evaluated at each grid-point
        Map_B: Numpy array; an NxM array containing the B-mode aperture mass evaluated at each grid-point
        Map_V: Numpy array; an NxM array containing the variance in the aperture mass evaluated at each grid-point
        kwargs; dict; keyword-arguments to be passed to photutil's detect_sources (default is SN>2.5; 25-pxs
    
    Returns:
        peak_catalog_E: Astropy Table; a table storing information about each peak in the E-mode
        peak_catalog_B: Astropy Table; a table storing information about each peak in the B-mode
    
    '''
    
    # first build the SN for the E/B-Modes
    sn_E = Map_E/np.sqrt(Map_V)
    sn_B = Map_B/np.sqrt(Map_V)
    
    # detect sources on each map
    seg_E = detect_sources(sn_E,**kwargs)
    seg_B = detect_sources(sn_B,**kwargs)
    
    
    # collecting interesting statistics from E-mode

    
    # collecting interesting statistics from E-mode, if any peaks exist
    if seg_E is not None:
    
        colnames = list(sn_columns.keys())
        output_colnames = list(sn_columns.values())
        catalog_sn_E = SourceCatalog(data=sn_E,segment_img=seg_E).to_table(colnames)
        peak_catalog_sn_E = catalog_sn_E.rename_columns(colnames,output_colnames)
        
        colnames = list(map_columns.keys())
        output_colnames = list(map_columns.values())
        catalog_ap_E = SourceCatalog(data=Map_E,segment_img=seg_E).to_table(colnames)
        peak_catalog_ap_E = catalog_ap_E.rename_columns(colnames,output_colnames)
        
        peak_catalog_E = hstack([catalog_sn_E,catalog_ap_E])
    
    else:
        
        # if there are no-sources in the B-mode, return an empty table
        peak_catalog_E = Table()

        
    # collecting interesting statistics from B-mode, if any peaks exist
    if seg_B is not None:
        
        colnames = list(sn_columns.keys())
        output_colnames = list(sn_columns.values())
        catalog_sn_B = SourceCatalog(data=sn_B,segment_img=seg_B).to_table(colnames)
        peak_catalog_sn_B = catalog_sn_B.rename_columns(colnames,output_colnames)
    
        colnames = list(map_columns.keys())
        output_colnames = list(map_columns.values())
        catalog_ap_B = SourceCatalog(data=Map_B,segment_img=seg_B).to_table(colnames)
        peak_catalog_ap_B = catalog_ap_B.rename_columns(colnames,output_colnames)
        
        peak_catalog_B = hstack([catalog_sn_B,catalog_ap_B])
        
    else:
        
        # if there are no-sources in the B-mode, return an empty table
        peak_catalog_B = Table()
    
    return peak_catalog_E, peak_catalog_B


def write_to_fits(Map_E,Map_B,Map_V,wcs,output_filename):
    '''
    This function write the mass-aperture statistics to a multi-extension fits
    
    Args:
        Map_E: Numpy array; an NxM array containing the E-mode aperture mass evaluated at each grid-point
        Map_B: Numpy array; an NxM array containing the B-mode aperture mass evaluated at each grid-point
        Map_V: Numpy array; an NxM array containing the variance in the aperture mass evaluated at each grid-point
        wcs; Astropy WCS; a wcs defined according to each grid-point
    
    Returns:
        None? Not sure what the outputs will look like yet

    '''
    
    # create the wcs header, pass it to the primary HDU
    header = wcs.to_header()
    
    Map_E_hdu = fits.PrimaryHDU(data=Map_E,header=header)
    Map_B_hdu = fits.ImageHDU(data=Map_B,name='M_B')
    Map_V_hdu = fits.ImageHDU(data=Map_V,name='M_V')
    
    # create the hdul and write it to disk
    hdul = fits.HDUList([Map_E_hdu,Map_B_hdu,Map_V_hdu])
    hdul.writeto(output_filename,overwrite=True)
    
    return


def draw_mass_map(Map_E,Map_B,Map_V,smoothing,Map_wcs):
    '''
    
    A helper function to render our nice mass_map plots
    
    Args:
        Map_E: MxN array; mass aperture map E-mode
        Map_B: MxN array; mass aperture map B-mode
        Map_V: MxN array; mass aperture map variance
        y_grid_samples: array; array containing the y-coordinates where Map is evaluated
        x_grid_samples: array; array containing the x-coordinates where Map is evaluated
        smoothing: float; the smoothing used for this mass map
    Returns:
        fig: matplotlib Figure; figure containing the plot
        (ax1,ax2): tuple of axes; a tuple containing the matplotlib axes
    '''
    
    fig,(ax1,ax2) = pl.subplots(1,2,figsize=(10,6),subplot_kw={'projection':Map_wcs,'coord_meta':{'unit':(u.deg,u.deg)},'aspect':'equal'})
    
    lon = ax1.coords[0]
    lat = ax1.coords[1]
    lon.set_major_formatter('d.d')
    lat.set_major_formatter('d.d')
    
    lon = ax2.coords[0]
    lat = ax2.coords[1]
    lon.set_major_formatter('d.d')
    lat.set_major_formatter('d.d')
    
    im = ax1.imshow(Map_E/np.sqrt(Map_V), vmin=-3, vmax=6,origin='lower')
    im = ax2.imshow(Map_B/np.sqrt(Map_V), vmin=-3, vmax=6,origin='lower')
    cbar = fig.colorbar(im, ax=(ax1,ax2), orientation="horizontal")
    cbar.ax.set_xlabel("S/N")
    ax1.set_xlabel('RA (deg)')
    ax1.set_ylabel('RA (deg)')
    ax2.set_xlabel('DEC (deg)')
    ax2.set_ylabel('DEC (deg)')
    
    ax1.set_title('E-Mode (Tangent)')
    ax2.set_title('B-Mode (Cross)')

    fig.suptitle(r"Aperture Mass S/N Map with $R_{\rm{ap}}$ = %.0f px"%(smoothing))
    
    return fig, (ax1,ax2)

    
# Hmm, the most correct way of doing this is as follows:
# 1. center mass_map on the X-ray peak/cluster-center
# 2. Define a LxL (L in degrees) on-sky
# 3. Compute distances w. a flat-sky; delta_theta = sqrt(delta_dec^2 + delta_ra^2*cos(dec)^2)
# 4. Define an appropriate WCS for this for saving!


if __name__ == '__main__':
    
    # collecting arguments from cln
    if len(sys.argv) == 10:

        table_filename = sys.argv[1]
        sample_spacing = float(sys.argv[2]) # leave in pixels
        Map_length = float(sys.argv[3]) # leave in degrees
        output_directory = sys.argv[4]
        cores = int(sys.argv[5])
        instrument = sys.argv[6]
        apply_cuts = int(sys.argv[7])
        cln = sys.argv[8]
        filter_name = sys.argv[9]
        
    else:
    
        print("python mass_map.py table coadd grid_resolution filter output_directory cores [OPTIONAL: LOWER_PATCH]")
        raise Exception("Improper Usage! Correct usage: python mass_map.py table coadd grid_resolution filter output_directory cores [OPTIONAL: LOWER_PATCH]")
    
    # first load the WCS to use while building Map
    butler = Butler('repo/repo')
    skymap = butler.get('skyMap',collections='skymaps',dataId={'instrument':'DECam','tract':0,'skymap':'{CLN}_skymap'.format(CLN=cln)})
    
    for tract in skymap:
        print(skymap.config)
    
    wcs = WCS(tract.wcs.getFitsMetadata())
    
    # read in the table, apply quality cuts, save the cut-table w/ an additional tag
    table = ascii.read(table_filename)
    
    if apply_cuts:
        table = load_quality_cuts(table,quality_cuts=quality_cuts)
        basename = Path(table_filename).stem
        table.write(output_directory + basename + '_Map_cut.csv',format='ascii.csv',overwrite=True)
    
    resolution = instrument_resolution[instrument]
    
    # load the aperture/filter
    q_filter = implemented_filters[filter_name]
    
    # load information about the cluster
    ned_result = Ned.query_object(cln)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    zL = ned_result[0]['Redshift']
    
    # define a grid centered on the cluster
    Map_center = skycoord_to_pixel(SkyCoord(ra=ra_cl,dec=dec_cl,unit='degree'),wcs=wcs)
    px_length = Map_length*3600/0.263
    x_sample = np.arange(Map_center[0] - px_length//2, Map_center[0] + px_length//2,sample_spacing)
    y_sample = np.arange(Map_center[1] - px_length//2, Map_center[1] + px_length//2,sample_spacing)
    y_grid,x_grid = np.meshgrid(y_sample,x_sample)
    
    # define a WCS for the sample-grid
    # originally built this block (w. Prakruth's help) for mass_maps for A360 w. ComCam
    Map_wcs = WCS(naxis=2)
    crval_sky = [ra_cl*u.deg,dec_cl*u.deg]
    Map_wcs.wcs.crval = [ra_cl,dec_cl]
    Map_wcs.wcs.crpix = [len(x_sample) // 2,len(y_sample) // 2] # center on the cluster
    delta_per_px = sample_spacing*resolution/3600
    Map_wcs.wcs.cdelt = [-delta_per_px,delta_per_px]
    Map_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    Map_wcs.wcs.radesys = 'ICRS'
    Map_wcs.wcs.equinox = 2000
    Map_wcs.wcs.cd = [[-delta_per_px,0],[0,delta_per_px]]
    
    # now that the grid is defined, collect information from the table
    x = table['x'].value
    y = table['y'].value
    g1 = table['g1'].value
    g2 = table['g2'].value
    
    # I've patched this, but I'll leave it here just in-case as a catch-all
    if 'shape_weight' not in table.colnames:
        weights = np.ones(len(x))/(0.365**2) # shape-weight is an inverse-variance, from previous studies e_rms ~ 0.365 
    else:
        weights = table['shape_weight'].value
    
    # temporary helper function for parallelizing over Rs
    def compute_and_measure_parallel(aperture_size):
        e_ap,b_ap,v_ap = compute_mass_map(x_grid,y_grid,x,y,g1,g2,weights,q_filter,filter_kwargs={'aperture_size':aperture_size})
        e_table,b_table = detect_mass_peaks(e_ap,b_ap,v_ap)
        return e_ap,b_ap,v_ap,e_table,b_table
    
    # distribute over apertures
    #TODO this can take a while and we can only really request 20-cores per-job that python can use, how can we distribute across nodes?
    with Pool(cores) as p:
        outputs = p.map(compute_and_measure_parallel,aperture_sizes)
    
    Map_E_table = Table()
    Map_B_table = Table()
    
    # unpack the outputs
    # I think unpacking in this manner is acceptable...
    # but I will wrap the drawing mass-map function into a separate call
    for i in range(len(aperture_sizes)):
        
        smoothing = aperture_sizes[i]
        
        Map_E = outputs[i][0]
        Map_B = outputs[i][1]
        Map_V = outputs[i][2]
        
        temp_E = outputs[i][3]
        temp_B = outputs[i][4]
        
        
        # Draw pretty pictures first
        fig, (ax1, ax2) = draw_mass_map(Map_E,Map_B,Map_V,smoothing,Map_wcs)
        
        # then update the table and draw the peaks
        # E-Mode
        if len(outputs[i][3]) != 0:
        
            e_table = outputs[i][3]
            e_table['aperture_size'] = np.ones(len(e_table)) * smoothing
            peak_skycoord = pixel_to_skycoord(xp=e_table['x_sn_max'],yp=e_table['y_sn_max'],wcs=Map_wcs)
            e_table['ra'] = peak_skycoord.ra.deg
            e_table['dec'] = peak_skycoord.dec.deg
            ax1.scatter(x=e_table['x_sn_max'], y=e_table['y_sn_max'],s=10,marker='x')
            e_table['SourceID'] = np.char.add(e_table['SourceID'].astype(str),'_' + str(smoothing)) # append the schirmer radius to make a unique source-ID
            Map_E_table = vstack([Map_E_table,e_table])
        
        #B-Mode
        if len(outputs[i][4]) != 0:
        
            b_table = outputs[i][4]
            b_table['aperture_size'] = np.ones(len(b_table)) * smoothing
            peak_skycoord = pixel_to_skycoord(xp=b_table['x_sn_max'],yp=b_table['y_sn_max'],wcs=Map_wcs)
            b_table['ra'] = peak_skycoord.ra.deg
            b_table['dec'] = peak_skycoord.dec.deg
            ax2.scatter(x=b_table['x_sn_max'],y=b_table['y_sn_max'],s=10,marker='x')
            b_table['SourceID'] = np.char.add(b_table['SourceID'].astype(str),'_' + str(smoothing)) # append the schirmer radius to make a unique source-ID
            Map_B_table = vstack([Map_B_table,b_table])
        
        # and save the outputs
        fig.savefig(output_directory + f'Map_SN_Rs_{smoothing}.png',bbox_inches='tight',dpi=720)
        write_to_fits(Map_E,Map_B,Map_V,Map_wcs,output_directory + f'Map_{filter_name}_{smoothing}.fits')
        
    # and lastly save the master Map tables
    Map_E_table.write(output_directory + f'Map_E_peak_catalog_{filter_name}.csv',format='ascii.csv',overwrite=True)
    Map_B_table.write(output_directory + f'Map_B_peak_catalog_{filter_name}.csv',format='ascii.csv',overwrite=True)
    
    
    
