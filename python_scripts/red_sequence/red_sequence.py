import sys
import os
import math
import treecorr

from pathlib import Path
from itertools import cycle

from multiprocessing import Pool

import numpy as np

from scipy.integrate import trapz, cumtrapz
from scipy.optimize import curve_fit, least_squares

from scipy import interpolate
from scipy import stats
from scipy.stats import norm, gaussian_kde, binned_statistic_2d

import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib.patches import Ellipse

import astropy.units as u

from astropy.io import ascii, fits
from astropy.table import Table, hstack, vstack
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip

from astroquery.ipac.ned import Ned

from lsst.daf.butler import Butler

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)

# homebrew modules below
from ..configs.red_sequence_config import instrument_resolution, filter_map, error_tag, color_bins_gr, color_bins_ri, mag_bins, zspec_tolerance, scale_upper_mag, z_pedestal, sn_cut
from ..misc.match_catalogs import match_catalogs

'''

This script slects members of the red-sequence by fitting to peaks in the CMD in g-r and r-i w.r.t. i-magnitude as a reference. Works remarkably well and borrows all the features used for LV II (bandwidth specified in kpc, variable upper-mag depending on redshift, etc).

'''

# two quick helper functions for running the linear fits
def model(p,x,y):
    return p[0] * np.array(x) + p[1] - np.array(y)
def model0(x,a,b):
    return a*x + b


# old implementation which does a rough-selection of the RS
def fit_RS_for_colors(gal_table,output_directory,color_0='g_cmodel_mag',color_1='r_cmodel_mag',mag='r_cmodel_mag',fitting_bins=np.linspace(15,18,15+1),mag_bounds=(15,18),mag_0_bound=22,member_tolerance=0.15,select_member_0=22,select_member_1=19):
    '''
    Right now our fitting is pretty simple; apply color cuts and fit a linear function to identifying the red-sequence.
    
    
    Args:
        gal_table: Astropy Table; a table of G galaxies to search for red-sequence members in
        output_directory: string; the filepath leading to an output directory for plots
        color_0: string; a string specifying the column in gal_table to use for the color (color_0 - color_1)
        color_1: string; a string specifying the column in gal_table to use for the color 
        mag: string; the reference magnitude to use in gal_table when searching for rs-members
        fitting_bins: Numpy array; an (N+1,) array specifying N magnitude-bins to search for the rs in
        mag_bounds: list-like; bounds on the magnitude to use for fitting the rs
        mag_0_bound: float; upper-bound for color_0
        member_tolerance: float; closeness to the RS best-fit for a galaxy to be considered a member
        select_member_0: float; maximum color_0 allowed for rs-members
        select_member_1: float; maximum color_1 allowed for rs-members
    
    Returns:
        rs_members: Astropy Table; a table of galaxies which are part of the RS
        select_rs; Numpy Array; a (G,)-length array which is true if the corresponding index in gal_table is in the rs
        
    '''
    
    # first apply magnitude cuts and prune NaN's
    select_fitting = (gal_table[mag] > mag_bounds[0]) & (gal_table[mag] < mag_bounds[1]) & (gal_table[color_0] < mag_0_bound)
    fit_table = gal_table[select_fitting]
    fit_table = fit_table[np.isfinite(fit_table[color_0]) & np.isfinite(fit_table[color_1]) & np.isfinite(fit_table[mag])]
    
    # setup a subplot for drawing the diagram
    fig,ax = pl.subplots()
    
    # and collect the relevant quantities for the fit
    colors = fit_table[color_0] - fit_table[color_1]
    ref_mag = fit_table[mag]
    
    # render the selected objects on a color-magnitude diagram
    ax.scatter(ref_mag,colors,marker='.',s=3,alpha=0.4)
    
    # first let's bin the objects in ref_mag
    median_color = []
    midpoints = (fitting_bins[1:] + fitting_bins[:-1])/2
    include_bin = []
    for i in range(len(fitting_bins) - 1):
        in_bin = (ref_mag > fitting_bins[i]) & (ref_mag < fitting_bins[i+1])
        clipped = sigma_clip(colors[in_bin],sigma=2,masked=False)
        median_color.append(np.median(clipped))
        include_bin.append(len(clipped) > 15) # if fewer than 15 objects exist in this particular slice/bin, drop it from the fitting
    
    median_color = np.array(median_color)
    
    # now run the fitting
    fit_median = median_color[include_bin]
    fit_mag = midpoints[include_bin]
    
    linear = lambda x, A, B : A*x + B
    popt, pcov = curve_fit(linear, fit_mag, fit_median)
    
    # draw the fitted curve on the color-magnitude diagram
    ax.plot(fit_mag,linear(fit_mag,*popt),'--',alpha=0.5)
    
    # select members
    select_rs = (gal_table[color_0] < select_member_0) & (gal_table[color_1] < select_member_1) 
    select_rs &= (np.abs((gal_table[color_0] - gal_table[color_1]) - linear(gal_table[mag],*popt)) < member_tolerance)
    rs_members = gal_table[select_rs]
    
    # highlight the members on the diagram
    ax.scatter(rs_members[mag],(rs_members[color_0] - rs_members[color_1]),marker='.',s=3)
    
    # lastly, make the plot pretty and save it
    ax.grid(ls=':')
    ax.set_xlim((14,19))
    ax.set_ylim((-1,3))
    ax.set_xlabel(mag[0]) # first-character is assumed the band
    ax.set_ylabel(color_0[0] + ' - ' + color_1[0])
    fig.savefig(output_directory + "%s%s_%s_diagram.png"%(color_0[0],color_1[0],mag[0]),bbox_inches='tight',dpi=720)
    
    return rs_members, select_rs


def fit_RS_color_binning(gal_table,output_filepath,filter_map=filter_map,color_0='g',color_1='r',ref_mag='i',mag_bins=np.linspace(13.5,19.5,21),mag_cut=20,color_bins=np.linspace(0.65,1,7),tolerance=0.09,min_per_bin=10,f_scale=0.09):
    '''
    A more sophisticated approach to identifying the red-sequence; first bin in magnitude then bin in color and select an overdensity in the color bin to use as a point for fitting the red-sequence
    
    
    Args:
        gal_table: Astropy Table; a table of G galaxies to search for red-sequence members in
        output_filepath: string; the filepath of an output quality check figure
        filter_map: dict; filter kde specifying column names given a band
        color_0: string; a string specifying the column in gal_table to use for the color (color_0 - color_1)
        color_1: string; a string specifying the column in gal_table to use for the color 
        ref_mag: string; the reference magnitude to use in gal_table when searching for rs-members
        mag_bins: Numpy array; an (N+1,) array specifying N magnitude-bins for identifying the RS
        mag_cut: float; the limit of RS galaxies, may vary per-cluster
        color_bins: Numpy array; an (M+1,) array specifying M color-bins for identifying the RS
        tolerance: float; nearness to the RS in color to be considered a member
        min_per_bin: int; minimum number of galaxies per-bin to select for fitting
        f_scale: float; scale used with the loss function
    
    Returns:
        rs_members: Astropy Table; a table of galaxies which are part of the RS
        rs_cut; Numpy Array; a (G,)-length array which is true if the corresponding index in gal_table is in the rs
        
    '''
    
    # let's collect the data needed from gal_table
    color = gal_table[filter_map[color_0]] - gal_table[filter_map[color_1]]
    mag = gal_table[filter_map[ref_mag]]
    
    # everything here can be done by using binned_statstic_2d routine in scipy
    count,mag_edges,color_edges,_ = binned_statistic_2d(x=mag,y=color,values=mag,statistic='count',bins=(mag_bins,color_bins))
    
    # first center the bins
    mag_bins = (mag_edges[1:] + mag_edges[:-1])/2
    color_bins = (color_edges[1:] + color_edges[:-1])/2
    
    # now compute the maximum array-wise down axis-1
    maximum_count = np.max(count,axis=1)
    maximum_arg = np.argmax(count,axis=1)
    
    #TODO there are some proper statistics we can do to select "overdensities" as +N-sigma variations
    # but that requreis knowing the distribution of galaxies w/out the cluster in the way
    # in principle we can kinda do this w. COSMOS and other datasets; but in-practice the mode should be robust enough
    # select the corresponding color bins
    select_bins = maximum_count > min_per_bin
    color_points = color_bins[maximum_arg[select_bins]]
    mag_points = mag_bins[select_bins]
    
    # alright, now that I have the points selected from the 2D binning we can run through 
    #TODO why this particular loss instead of doing a weighted LSQ using the per-bin scatters?
    result = least_squares(model, x0=np.array([0,1]),loss='soft_l1',f_scale=f_scale,args=(mag_points,color_points))
    
    # lastly, select red-sequence members using this model!
    # be sure to constrain it to the range of magnitudes 
    expected_color = model0(mag, *result.x)
    rs_cut = np.abs(expected_color - color) < f_scale
    rs_cut &= (mag < mag_cut) & (mag > mag_edges[0])
    
    # alright, right now I have all the bins collected together along with RS-members selected, let's draw a quality-check plot!
    # two-panel, one with the CMD and another the histogram per mag-bin.
    fig,(ax1,ax2) = pl.subplots(1,2,figsize=(6.4*2,4.8))
    ax1.scatter(x=mag,y=color,s=1,label='Galaxy')
    ax1.scatter(x=mag[rs_cut],y=color[rs_cut],s=20,color='C3',alpha=0.2,label='Candidate RS-Member')
    ax1.scatter(x=mag_points,y=color_points,s=25,color='C2',label='Binned Peaks')
    ax1.plot(mag_bins,model0(mag_bins, *result.x),color='C2',linestyle='--')
    ax1.set_ylim((color_bins[0]/1.25,color_bins[-1]*1.25))
    ax1.set_xlim((mag_bins[0]/1.05,mag_bins[-1]*1.05))
    ax1.legend(loc='upper left')
    # put together everything for a colormap on ax2
    cmap = mpl.cm.viridis_r
    norm = mpl.colors.Normalize(vmin=19.5,vmax=13.5)
    colors = pl.cm.viridis((mag_bins - np.min(mag_points))/np.max(mag_bins - np.min(mag_points)))
    
    for i in range(len(mag_points)):
        if select_bins[i]:
            draw_me = count[i,:]
            draw_me = draw_me/np.sum(draw_me)
            ax2.plot(color_bins,draw_me,alpha=0.7,color=colors[i])
    
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax2)
    cb.set_label(f'{ref_mag}',rotation='horizontal')
    ax1.set_xlabel(f'{ref_mag}')
    ax1.set_ylabel(f'{color_0} - {color_1}')
    ax2.set_xlabel(f'{color_0} - {color_1}')
    ax2.set_ylabel('Frequency')
    slope = result.x[0]
    intercept = result.x[1]
    text = f'$({color_0} - {color_1}) = m~({ref_mag}) + b $ \n $ m = {slope:.3f}$ \n $ b = {intercept:1.3f} $'
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax1.text(x=0.05,y=0.05,s=text,transform=ax1.transAxes,bbox=props,ma='center')
    pl.tight_layout()
    
    # finally, plotting is done and I can save this figure
    fig.savefig(output_filepath,dpi=720,bbox_inches='tight')
    
    # and I can return the outputs
    rs_members = gal_table[rs_cut]
    return rs_members, rs_cut


#TODO if we add more kernels, we should move them to a separate script
def gaussian_circular_2d(x,y,x0=0,y0=0,sigma=1):
    return (1/(2*np.pi*(sigma**2))) * np.exp( -0.5*( (x - x0)**2 + (y-y0)**2 )/sigma**2 )


#TODO like w. mass_map, there is probably a clever way of avoiding the nested for-loop here w. broadcasting
#TODO also again eventually we should move this to an on-sky method
def build_gaussian_kde(rs_member_x,rs_member_y,x_grid,y_grid,weights=None,kernel=gaussian_circular_2d,kernel_kwargs={'sigma':100}):
    '''
    This function evaluates a KDE sampled on a given grid.
    
    Args:
        rs_member_x: Numpy array; an (n,)-shaped array of x-coordinates for each object
        rs_member_y: Numpy array; an (n,)-shaped array of y-coordinates for each object
        x_grid: Numpy array; an NxM array of x-coordinates to sample the aperture-mass on
        y_grid: Numpy array; an NxM array of y-coordinates to sample the aperture-mass on
        weights: Numpy array; an (n,)-shaped array of per-object weights. If None, defaults to uniform weight
        kernel: a kernel to use for building the KDE
        kernel_kwargs: a dictionary of kwargs to pass to kernel
    
    Returns:
        kde: Numpy Array; an NxM array containing the KDE at each point
    
    '''
    
    if weights is None:
        weights = np.ones(len(rs_member_x))
    
    y_shape = len(y_grid[:,0])
    x_shape = len(x_grid[0,:])
    
    # create an array to store the kde
    kde = np.zeros((y_shape,x_shape))
    
    for i in range(y_shape):
        for j in range(x_shape):
            
            kde_values = kernel(x=x_grid[j,i],y=y_grid[j,i],x0=rs_member_x,y0=rs_member_y,**kernel_kwargs)
            kde[i,j] = np.sum(kde_values*weights)
    
    return kde

def write_to_fits(kde,wcs,output_filename):
    '''
    This function write the mass-aperture statistics to a multi-extension fits
    
    Args:
        kde: Numpy array; an NxM array containing the kde of the number density evaluated at each grid-point
        wcs; Astropy WCS; a wcs defined according to each grid-point
    
    Returns:
        None? Not sure what the outputs will look like yet

    '''
    
    # create the wcs header, pass it to the primary HDU
    header = wcs.to_header()
    
    kde_hdu = fits.PrimaryHDU(data=kde,header=header)
    
    # create the hdul and write it to disk
    hdul = fits.HDUList([kde_hdu])
    hdul.writeto(output_filename,overwrite=True)
    
    return

if __name__ == '__main__':
    
    # collecting arguments from cln
    if len(sys.argv) == 8:
    
        gal_table_filename = sys.argv[1]
        sample_spacing = int(sys.argv[2]) # leave in px
        kde_length = float(sys.argv[3]) # leave in degrees
        bandwidth = int(sys.argv[4]) # in kpc
        output_directory = sys.argv[5]
        instrument = sys.argv[6]
        cluster_name = sys.argv[7]
        specz_cat_filename = None
        specz_header = None
    
    elif len(sys.argv) == 10:
        
        gal_table_filename = sys.argv[1]
        sample_spacing = int(sys.argv[2]) 
        kde_length = float(sys.argv[3]) # leave in degrees
        bandwidth = int(sys.argv[4]) # in kpc
        output_directory = sys.argv[5]
        instrument = sys.argv[6]
        cluster_name = sys.argv[7]
        specz_cat_filename = sys.argv[8]
        specz_header = sys.argv[9] 
    
    else:
    
        raise Exception("Improper Usage! Correct usage: python red_sequence.py gal_table coadd grid_resolution bandwidth output_directory instrument lower_patch")
    
    # load the wcs from the skyMap
    butler = Butler('repo/repo')
    skyMap = butler.get('skyMap',collections='skymaps',dataId={'instrument':'DECam','tract':0,'skymap':'{CLN}_skymap'.format(CLN=cluster_name)})
    
    for tract in skyMap:
        print(skyMap.config)
    
    wcs = WCS(tract.wcs.getFitsMetadata())
    
    # collect the cluster info
    ned_result = Ned.query_object(cluster_name)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    zL = ned_result[0]['Redshift']
    
    # load the instrument resolution
    resolution = instrument_resolution[instrument]
    
    # if enabled, scale the upper magnitude based on the pedital (z~0.08) by default
    if scale_upper_mag:
        min_mag = mag_bins[0]
        max_mag = mag_bins[-1]
        distance_ratio = cosmo.luminosity_distance(zL)/cosmo.luminosity_distance(z_pedestal)
        max_mag = max_mag + 5*np.log10( distance_ratio.value )
        mag_bins = np.linspace(min_mag,max_mag,len(mag_bins))
    
    gal_table = ascii.read(gal_table_filename)
    
    # get ready to filter out low-sn detections
    cut_sn = (1/gal_table[filter_map['r'] + error_tag] > sn_cut*np.log(10)/2.5)
    
    # just in-case, filter any NaN's in g/r/i
    cut_nan = np.isfinite(gal_table[filter_map['r']]) & np.isfinite(gal_table[filter_map['g']]) & np.isfinite(gal_table[filter_map['i']])
    cut_nan &= np.isfinite(gal_table[filter_map['r'] + error_tag]) & np.isfinite(gal_table[filter_map['g'] + error_tag]) & np.isfinite(gal_table[filter_map['i'] + error_tag])
    
    # and lastly select w/in 1-deg of the cluster center
    coord = SkyCoord(ra=gal_table['ra'],dec=gal_table['dec'],unit='degree')
    cut_center = coord.separation(SkyCoord(ra=ra_cl,dec=dec_cl,unit='degree')) < 60*u.arcmin
    
    # actually apply the cuts
    cut_gal_table = gal_table[cut_sn & cut_nan & cut_center]
    
    # define a grid centered on the cluster
    kde_center = skycoord_to_pixel(SkyCoord(ra=ra_cl,dec=dec_cl,unit='degree'),wcs=wcs)
    px_length = kde_length*3600/0.263
    x_sample = np.arange(kde_center[0] - px_length//2, kde_center[0] + px_length//2,sample_spacing)
    y_sample = np.arange(kde_center[1] - px_length//2, kde_center[1] + px_length//2,sample_spacing)
    y_grid,x_grid = np.meshgrid(y_sample,x_sample)
    
    # define a WCS for the sample-grid
    # originally built this block (w. Prakruth's help) for mass_kdes for A360 w. ComCam
    kde_wcs = WCS(naxis=2)
    crval_sky = [ra_cl*u.deg,dec_cl*u.deg]
    kde_wcs.wcs.crval = [ra_cl,dec_cl]
    kde_wcs.wcs.crpix = [len(x_sample) // 2,len(y_sample) // 2] # center on the cluster
    delta_per_px = sample_spacing*resolution/3600
    kde_wcs.wcs.cdelt = [-delta_per_px,delta_per_px]
    kde_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    kde_wcs.wcs.radesys = 'ICRS'
    kde_wcs.wcs.equinox = 2000
    kde_wcs.wcs.cd = [[-delta_per_px,0],[0,delta_per_px]]
    
    # use the peak-finding code for this, first initial selection in gr then final in ri
    gr_members,gr_selections = fit_RS_color_binning(cut_gal_table,output_directory + 'photom_gr_rs_selection.png',filter_map=filter_map,color_0='g',color_1='r',color_bins=color_bins_gr,tolerance=0.09,min_per_bin=15,mag_bins=mag_bins,f_scale=0.09)
    ri_members,ri_selections = fit_RS_color_binning(gr_members,output_directory + 'photom_ri_rs_selection.png',filter_map=filter_map,color_0='r',color_1='i',color_bins=color_bins_ri,tolerance=0.09,min_per_bin=15,mag_bins=mag_bins,f_scale=0.09)
    
    # if a spec-z catalog is provided, run this same analysis on a spec-z matched catalog containing just cluster members (truth) to check for purity and completeness
    # using while here to make it easy to leave this block,
    while specz_cat_filename is not None:
            
        # first load in the specz-table
        gal_table_specz = ascii.read(specz_cat_filename)
        
        # now filter for SN, NaN
        cut_sn = (1/gal_table_specz[filter_map['r'] + error_tag] > sn_cut*np.log(10)/2.5)
        cut_nan = np.isfinite(gal_table_specz[filter_map['r']]) & np.isfinite(gal_table_specz[filter_map['g']]) & np.isfinite(gal_table_specz[filter_map['i']])
        cut_nan &= np.isfinite(gal_table_specz[filter_map['r'] + error_tag]) & np.isfinite(gal_table_specz[filter_map['g'] + error_tag]) & np.isfinite(gal_table_specz[filter_map['i'] + error_tag])
        
        # and lastly select w/in 1-deg of the cluster center
        coord = SkyCoord(ra=gal_table_specz['ra'],dec=gal_table_specz['dec'],unit='degree')
        cut_center = coord.separation(SkyCoord(ra=ra_cl,dec=dec_cl,unit='degree')) < 60*u.arcmin
        
        # apply the cuts and check that there are enough objects
        cut_gal_table_specz = gal_table_specz[cut_sn & cut_nan & cut_center]
        if len(gal_table_specz) < 1000:
            print('Spec-z catalog exists, but has too few entries near the cluster for a reliable purity/completeness check!')
            break
        
        cluster_members = np.abs(cut_gal_table_specz[specz_header] - zL) < zspec_tolerance
        cut_gal_table_specz_members = cut_gal_table_specz[cluster_members]
        
        # run the RS- on selected cluster members
        try:
            gr_members_spec,gr_selections_spec = fit_RS_color_binning(cut_gal_table_specz_members,output_directory + 'spec_members_gr_rs_selection.png',filter_map=filter_map,color_0='g',color_1='r',color_bins=color_bins_gr,tolerance=0.09,min_per_bin=10,mag_bins=mag_bins,f_scale=0.09)
            ri_members_spec,ri_selections_spec = fit_RS_color_binning(gr_members_spec,output_directory + 'spec_members_ri_rs_selection.png',filter_map=filter_map,color_0='r',color_1='i',color_bins=color_bins_ri,tolerance=0.09,min_per_bin=10,mag_bins=mag_bins,f_scale=0.09)
            
            # to check for completeness/purity, run RS-selection on all entries in specz-catalog
            gr_members_spec_all,gr_selections_spec = fit_RS_color_binning(cut_gal_table_specz,output_directory + 'spec_all_gr_rs_selection_specz.png',filter_map=filter_map,color_0='g',color_1='r',color_bins=color_bins_gr,tolerance=0.09,min_per_bin=10,mag_bins=mag_bins,f_scale=0.09)
            ri_members_spec_all,ri_selections_spec = fit_RS_color_binning(gr_members_spec_all,output_directory + 'spec_all_ri_rs_selection_all_specz.png',filter_map=filter_map,color_0='r',color_1='i',color_bins=color_bins_ri,tolerance=0.09,min_per_bin=10,mag_bins=mag_bins,f_scale=0.09)
        
        except:
            print("Not enough confirmed cluster members to check for completeness/purity!")
            print("Skipping these statistics and rendering the RS-map")
            break
        
        #TODO when we have per-object identifiers, these objects can be matched by their ID rather than by an explicit coordinate matching
        # now that I have everythign collected together, its time to check the purity of these catalogs, first lets run a quick matching
        gr_matches = match_catalogs(gr_members_spec_all,'decam',gr_members_spec,'decam',refcat_tag='_photom1',catalog_tag='_spec2',sep=0.1*u.arcsec)
        ri_matches = match_catalogs(ri_members_spec_all,'decam',ri_members_spec,'decam',refcat_tag='_photom1',catalog_tag='_spec2',sep=0.1*u.arcsec)
        
        # now I can compute the TP/FP/FN number just from the lengths of these arrays
        tp_gr = len(gr_matches)
        tp_ri = len(ri_matches)
        fp_gr = len(gr_members_spec_all) - len(gr_matches) # false positives in gr_members_spec_all, but not in matched
        fp_ri = len(ri_members_spec_all) - len(ri_matches)
        fn_gr = len(gr_members_spec) - len(gr_matches) # false negatives are in gr_members_spec, but not in matched
        fn_ri = len(ri_members_spec) - len(ri_matches)
        
        gr_purity = tp_gr / (tp_gr + fp_gr) # of the "positives", what fraction are correctly assigned as positive?
        gr_completeness = tp_gr / (tp_gr + fn_gr) # what fraction of the "positives" were selected out of all "true positives"?
        ri_purity = tp_ri / (tp_ri + fp_ri)
        ri_completeness = tp_ri / (tp_ri + fn_ri) 
        
        # finally, lets draw a plot to summarize everything!
        fig,(ax1,ax2) = pl.subplots(1,2,figsize=(6.4*2,4.8))
        
        # first, for ax1 lets draw the gr stuff
        gr_all = cut_gal_table_specz[filter_map['g']] - cut_gal_table_specz[filter_map['r']]
        mag_all = cut_gal_table_specz[filter_map['i']]
        ax1.scatter(x=mag_all,y=gr_all,s=1,label='Background Galaxy')
        
        gr_all_rs = gr_members_spec_all[filter_map['g']] - gr_members_spec_all[filter_map['r']]
        mag_all_rs = gr_members_spec_all[filter_map['i']]
        ax1.scatter(x=mag_all_rs,y=gr_all_rs,s=15,color='C2',alpha=0.3,label='Candidate RS-Member')
        
        gr_true_rs = gr_members_spec[filter_map['g']] - gr_members_spec[filter_map['r']]
        mag_true_rs = gr_members_spec[filter_map['i']]
        ax1.scatter(x=mag_true_rs,y=gr_true_rs,s=25,color='C3',alpha=0.5,label='True RS-Member',marker='x')
        
        # now setup the labels/legends for ax1, fortunately I have a blueprint for this already
        ax1.legend(loc = 'upper left')
        ax1.set_xlabel('i')
        ax1.set_ylabel('g - r')
        ax1.set_xlim(( mag_bins[0]/1.05, mag_bins[-1]*1.05 ))
        ax1.set_ylim(( color_bins_gr[0]/1.25, color_bins_gr[-1]*1.25 ))
        text = f' Completeness: {gr_completeness:.2f} \n Purity: {gr_purity:.2f} '
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax1.text(x=0.05,y=0.05,s=text,transform=ax1.transAxes,bbox=props,ma='center')
        
        # second, for ax2 lets draw the ri stuff
        ri_all = cut_gal_table_specz[filter_map['r']] - cut_gal_table_specz[filter_map['i']]
        mag_all = cut_gal_table_specz[filter_map['i']]
        ax2.scatter(x=mag_all,y=ri_all,s=1,label='Galaxy')
        
        ri_all_rs = ri_members_spec_all[filter_map['r']] - ri_members_spec_all[filter_map['i']]
        mag_all_rs = ri_members_spec_all[filter_map['i']]
        ax2.scatter(x=mag_all_rs,y=ri_all_rs,s=15,color='C2',alpha=0.3,label='Candidate RS-Member')
        
        ri_true_rs = ri_members_spec[filter_map['r']] - ri_members_spec[filter_map['i']]
        mag_true_rs = ri_members_spec[filter_map['i']]
        ax2.scatter(x=mag_true_rs,y=ri_true_rs,s=25,color='C3',alpha=0.5,label='True RS-Member',marker='x')
        
        # now setup the labels/legends for ax1, fortunately I have a blueprint for this already
        ax2.legend(loc='upper left')
        ax2.set_xlabel('i')
        ax2.set_ylabel('r - i')
        ax2.set_xlim(( mag_bins[0]/1.05, mag_bins[-1]*1.05 ))
        ax2.set_ylim(( color_bins_ri[0]/1.25, color_bins_ri[-1]*1.25 ))
        text = f' Completeness: {ri_completeness:.2f} \n Purity: {ri_purity:.2f} '
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax2.text(x=0.05,y=0.05,s=text,transform=ax2.transAxes,bbox=props,ma='center')
        
        # finally get a tight_layout running ands save
        pl.tight_layout()
        fig.savefig(output_directory + 'specz_rs_selection.png',dpi=720,bbox_inches='tight')
        
        # and leave the block
        break
        
    # alright, that's the quality check done! Next up, time to create the KDE from these
    rs_members = ri_members
    rs_member_x = rs_members['x']
    rs_member_y = rs_members['y']
    
    #TODO what is the best choice of weights for this?
    # bandwidth is in kpc at the cluster, convert this to a proper bandwidth in pixel-coordinates
    bandwidth = ( bandwidth / cosmo.angular_diameter_distance(zL).to_value('kpc') ) * ( 180/np.pi ) * 3600 / resolution
    
    kde = build_gaussian_kde(rs_member_x,rs_member_y,x_grid,y_grid,kernel_kwargs={'sigma':bandwidth})
    
    # upscale kde to reflect the sampling, e.g. kde is saved in density/spx as before
    kde = kde*(sample_spacing**2)

    write_to_fits(kde,wcs=kde_wcs,output_filename=output_directory + 'rs_density_kde.fits')
    
    # one-last thing, render the kde for reference
    fig,ax = pl.subplots(subplot_kw={'projection':kde_wcs,'coord_meta':{'unit':(u.deg,u.deg)},'aspect':'equal'})
    im = ax.imshow(kde,origin='lower')
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_major_formatter('d.d')
    lat.set_major_formatter('d.d')
    
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('DEC (deg)')
    
    cbar = fig.colorbar(im,ax=ax,orientation="vertical")
    cbar.ax.set_ylabel("Density-per-spx")
    fig.savefig(output_directory + 'rs_density.png',bbox_inches='tight',dpi=720)






