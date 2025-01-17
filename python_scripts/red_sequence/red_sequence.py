import sys
import os
import math
import treecorr

from pathlib import Path
from itertools import cycle

from multiprocessing import Pool

import numpy as np

from scipy.integrate import trapz, cumtrapz
from scipy.optimize import curve_fit

from scipy import interpolate
from scipy import stats
from scipy.stats import norm, gaussian_kde

import matplotlib.pyplot as pl
from matplotlib.patches import Ellipse

from astropy.io import ascii, fits
from astropy.table import Table, hstack, vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip

from astroquery.ipac.ned import Ned

#TODO overhaul this to a standard-IO
def load_resolution():
    '''
    A temporary function for storing a hand-coded dictionary specifying the resolution of different instruments
    
    Args:
        None
    
    Returns:
        instrument_reslution: Dictionary; a dictionary containing the resolution for keyed instruments
    
    '''
    
    instrument_resolution = {
                             'decam' : 0.263,
                             'hsc' : 0.168,
                            }
    
    return instrument_resolution

'''

This script runs a ROUGH measurement of the red-sequence galaxies and estimates the number-density of RS-members. The goal is to fit multiple linear-functions to the red-sequence in a color-magnitude diagram, then select member-galaxies which are sufficiently close to the curve. By default, we do this in two diagrams (g-r v. r & r-i v. i). Our latest analysis (LVII) does a somewhat more thorough analysis...
  
  1. For fitting the curves, select r(i) in (15,18) with g(r) < 22; use cmodel-mags and objects w. S/N > 10
  
  2. Draw two color-magnitude diagrams.
  
  3. Create magnitude-bins, assign median-color w. 2sig-clipping
  
  4. Include bins w. >15-objects in a linear-fit [ color(mag) ]; this is (roughly) the RS for a cluster
  
  5. Select objects w. | color - RS_color | < 0.15; these are roughly members.
  
  6. Generate and KDE+WCS for the density of objects

'''

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
    #TODO there is probably a nicer numpy-array like way of doing this
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
    if len(sys.argv) == 7:

        gal_table_filename = sys.argv[1]
        coadd_filename = sys.argv[2]
        sample_spacing = int(sys.argv[3])
        bandwidth = int(sys.argv[4])
        output_directory = sys.argv[5]
        instrument = sys.argv[6]
        lower_patch = (3,3)
        
    elif len(sys.argv) == 8:
    
        gal_table_filename = sys.argv[1]
        coadd_filename = sys.argv[2]
        sample_spacing = int(sys.argv[3]) 
        bandwidth = int(sys.argv[4])
        output_directory = sys.argv[5]
        instrument = sys.argv[6]
        lower_patch = sys.argv[7].split(',')
    
    else:
    
        raise Exception("Improper Usage! Correct usage: python red_sequence.py gal_table coadd grid_resolution bandwidth output_directory instrument [OPTIONAL: LOWER_PATCH]")
    
    sn_cut = 10
    
    gal_table = ascii.read(gal_table_filename)
    cut_sn = (1/gal_table['r_cmodel_magerr'] > sn_cut*np.log(10)/2.5)
    cut_gal_table = gal_table[cut_sn]
    
    # load the wcs from the coadd
    header = fits.getheader(coadd_filename,0)
    coadd_wcs = WCS(header)
    coadd_shape = coadd_wcs.array_shape
    
    # kde_shape should be (coadd[0]/res,coadd[1]/res)
    kde_shape = (int(coadd_shape[0]/sample_spacing),int(coadd_shape[1]/sample_spacing))
    
    # build a grid centered on the coadd
    centering_y = coadd_shape[0] - (kde_shape[0]*sample_spacing)
    centering_x = coadd_shape[1] - (kde_shape[1]*sample_spacing)
    
    #TODO what is the best way of telling kde the patch this coadd is located in (which offsets the px-coordinates by 4k px/patch). I've added an option to pass this via cln but, there may be a better solution?
    #TODO we could add information on the patch to the header of the coadd, that way we don't need to call the lower-patch nor the px/patch explicitly here nor at lines further below
    x_grid_samples = np.arange(int(lower_patch[1])*4000, (int(lower_patch[1])*4000) + coadd_shape[1],sample_spacing) + sample_spacing/2 + centering_x/2
    y_grid_samples = np.arange(int(lower_patch[0])*4000, (int(lower_patch[0])*4000) + coadd_shape[0],sample_spacing) + sample_spacing/2 + centering_y/2
    y_grid,x_grid = np.meshgrid(y_grid_samples,x_grid_samples)
    
    gr_selected,gr_cut = fit_RS_for_colors(cut_gal_table,output_directory)
    ri_selected,ri_cut = fit_RS_for_colors(cut_gal_table,output_directory,color_0='r_cmodel_mag',color_1='i_cmodel_mag',mag='i_cmodel_mag')
    
    rs_members = cut_gal_table[gr_cut & ri_cut]
    rs_member_x = rs_members['x']
    rs_member_y = rs_members['y']
    
    kde = build_gaussian_kde(rs_member_x,rs_member_y,x_grid,y_grid,kernel_kwargs={'sigma':bandwidth})
    
    # upscale kde to reflect the sampling, e.g. kde is saved in density/spx as before
    kde = kde*(sample_spacing**2)
    
    #TODO same problem from mass_map, these should probably be cleaned-up and wrapped in a function
    resolution_dict = load_resolution()
    resolution = resolution_dict[instrument] # ("/px)
    
    # create the wcs for kde; subtract the lower-patch-index since the coadd-wcs has px 4000*LOWER_INDEX => 0
    # use pixel 0,0 as the crval
    kde_wcs = WCS(naxis=2)
    crval_sky = coadd_wcs.pixel_to_world(x_grid_samples[0]-int(lower_patch[1])*4000,y_grid_samples[0]-int(lower_patch[0])*4000)
    kde_wcs.wcs.crval = [crval_sky.ra.degree,crval_sky.dec.degree]
    kde_wcs.wcs.crpix = [0,0]
    kde_wcs.wcs.cdelt = [-sample_spacing*resolution/3600,sample_spacing*resolution/3600]
    kde_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    kde_wcs.wcs.radesys = 'ICRS'
    kde_wcs.wcs.equinox = 2000

    write_to_fits(kde,wcs=kde_wcs,output_filename=output_directory + 'rs_density_kde.fits')
    
    # one-last thing, render the kde for reference
    fig,ax = pl.subplots()
    im = ax.imshow(kde,origin='lower',extent=[np.min(x_grid_samples), np.max(x_grid_samples), np.min(y_grid_samples), np.max(y_grid_samples) ])
    ax.set_xlabel('x [px]')
    ax.set_ylabel('y [px]')
    
    cbar = fig.colorbar(im,ax=ax,orientation="vertical")
    cbar.ax.set_ylabel("Density-per-spx")
    fig.savefig(output_directory + 'rs_density.png',bbox_inches='tight',dpi=720)






