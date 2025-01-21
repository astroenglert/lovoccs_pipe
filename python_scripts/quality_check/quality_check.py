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

from astroquery.ipac.ned import Ned

# homebrew modules below
from ..configs.quality_check_config import filter_map, instrument_resolution, quality_cuts, possible_bands

'''
This is for checking the quality of our measurements (e.g. noise from atm/optics)... to do this we pass magnitudes which are:
  - zp-corrected (consistent w/ refcats)
  - not de-reddened (this isn't about photons absorbed while passing through the mily-way)
  - one star-catalog which is matched with gaia

There is no real algorithm for this... instead it's just about drawing a few pretty-pictures, this includes:

1. plot_mag_hist; model high-mag count with a power law, then extrapolate to where the histogram-count is half the power-law count
  - effectively, compute the magnitude where we hit 50% completeness
  
2. plot_snr_mag_hist; draw histograms of objects w. SN 9-11 and SN 19-21

3. plot_mag_magerr; scatter plot of magerr (or-SN directly?) v. mag

4. compare_astrometry; check the astrometry against Gaia

5. star_gal_correlation; check that the positions of stars/gals are uncorrelated,
  - use a method called treecorr
  - check correlation between e1/e2_star & e1/e2_gal (C1); e1_star & e2_gal (C2); e2_star & e1_gal (C3)

6. draw_psf; scatter plots of psf-stars + fwhm/e1/e2 histograms, reminscent of check_visit
  - I don't think fitting is appropriate here (since really at each px the coadd has a stacked-psf from visits), so kron-radius (first-radial-moment) for now... eventually we should use the catalogs from LSP to plot a better statistic
  
7. draw_visits; draw the visits covering the cluster
  - load from query_result_CLN

8. draw_catalog; draw the complete catalog

'''


#TODO overhaul this to a standard-IO or function for running quality-cuts read from a config (this is cp'ed from mass_fit for now...)
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
                        #('g1','<',2),
                        #('g2','<',2), # and these two together are equivalent to |g| < 2
                        #('chi2_mod','<',4), # beware, in our old catalogs this is named chi2_mod, new is mod_chi2
                        #('odds','>',0.95),
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
    
    return cut_table


def plot_mag_hist(star_catalog,gal_catalog,output_directory,fit_range=[23,24],mag_bins=np.linspace(22,27,101),skymap_center=(24000,24000),count_range=(5e2,3e4),max_sep=13500):
    '''
    Run a completeness test: Model the tail of n(m), then find where the histogram drops to below half that predicted by the tail-model.
    
    Args:
        star_catalog: Astropy Table; an Astropy table containing a ZP-corrected catalog of stars
        gal_catalog: Astropy Table; an Astropy table containing a ZP-corrected catalog of stars
        output_directory: string; a string specifying an output directory
        fit_range: array-like: the range of magnitudes used for modelling the tail
        mag_bins: Numpy array: an (N+1,) array specifying the edges of N magnitude bins
        skymap_center: list-like; a list specifying the px (x,y)-coordinates to use for a distance-cut w. max_sep
        count_range: list-like; a list specifying the y-lim of the output plot
        max_sep: float; the maximum separation for a source to be included in these plots
      
    Returns:
        None?
    
    '''
    
    col_names = star_catalog.colnames
    
    # lazy-way of getting the number of bands which exist in the catalog
    num_bands = 0
    bands_in_catalog = []
    for i in possible_bands: 
        if (i + '_psf_mag') in col_names:
            num_bands += 1
            bands_in_catalog.append(i)
    
    # defining a temporary linear function for fitting
    f_linear = lambda x, A, B: A*x + B
    
    # cut to the center ~1deg
    #TODO center of the skymap isn't always aligned with the center of the cluster/our-pointing (although it usually is...), is there a way we can fix this even for targets where we have a non-standard pointing?
    star_center_cut = (np.sqrt((star_catalog['x'] - skymap_center[0])**2 + (star_catalog['y']-skymap_center[1])**2) < max_sep)
    gal_center_cut = (np.sqrt((gal_catalog['x'] - skymap_center[0])**2 + (gal_catalog['y']-skymap_center[1])**2) < max_sep)
    
    star_center = star_catalog[star_center_cut]
    gal_center = gal_catalog[gal_center_cut]
    
    #TODO if we tweak the star/gal separation and have failed-classifications, we might drop objects we'd like to otherwise include by re-assembling the all-catalog
    all_center = vstack((star_center,gal_center))
    
    
    # now create the figure and start plotting!
    rows = math.ceil((num_bands)/2)
    fig, axs = pl.subplots(2, rows, figsize=(rows*3.7,7))
    axs[-1,-1].axis('off') # clear a space for the legend later!
    
    for i,band in enumerate(bands_in_catalog):
        
        # load the correct axes
        ax = axs[i//rows, i%rows]
        
        # collect the cmodel/psf-mags
        star_mag = all_center[band + '_psf_mag']
        gal_mag = all_center[band + '_cmodel_mag']
        star_mag = star_mag[np.isfinite(star_mag)]
        gal_mag = gal_mag[np.isfinite(gal_mag)]
        
        # build histograms
        hist_star, bin_edges_star = np.histogram(star_mag,bins=mag_bins)
        midpoint_star = (bin_edges_star[1:] + bin_edges_star[:-1])/2
        hist_gal, bin_edges_gal = np.histogram(gal_mag,bins=mag_bins)
        midpoint_gal = (bin_edges_gal[1:] + bin_edges_gal[:-1])/2
        
        # select range and fit the tail of the distribution
        select_star = (midpoint_star > fit_range[0]) & (midpoint_star < fit_range[1])
        select_gal = (midpoint_gal > fit_range[0]) & (midpoint_gal < fit_range[1])
        popt_star, pcov_star = curve_fit(f_linear, midpoint_star[select_star], np.log10(hist_star[select_star]))
        popt_gal, pcov_gal = curve_fit(f_linear, midpoint_gal[select_gal], np.log10(hist_gal[select_gal]))
        
        # the interpolation here helps smooth-out the limiting-mag
        interp_star = interpolate.interp1d(midpoint_star, np.log10(2*hist_star))
        interp_gal = interpolate.interp1d(midpoint_gal, np.log10(2*hist_gal))
        
        # not sure why this is called f_root... but I'm inheriting this naming from SF's script
        root_star = lambda x : interp_star(x) - f_linear(x,*popt_star)
        root_gal = lambda x : interp_gal(x) - f_linear(x,*popt_gal)
        
        # search for when the completeness drops to 50%
        completeness_mag_star = midpoint_star[np.argmin(np.abs(root_star(midpoint_star)))]
        completeness_mag_gal = midpoint_gal[np.argmin(np.abs(root_gal(midpoint_gal)))]
        
        # drawing the histograms and a vline at 50% completeness
        ax.plot(midpoint_star, hist_star, label="PSF Mag")
        ax.plot(midpoint_gal, hist_gal, label="CModel Mag")
        ax.axvline(x=completeness_mag_star, linestyle='--', color='C0', alpha=0.7)
        ax.axvline(x=completeness_mag_gal, linestyle='--', color='C1', alpha=0.7)
        
        # setting title, scale, etc
        ax.set_yscale('log')
        ax.set_title("%s: PSF %.1f, CModel %.1f"%(band, completeness_mag_star, completeness_mag_gal) )
        ax.set_xlim((np.min(mag_bins),np.max(mag_bins)))
        ax.set_ylim((count_range))
        ax.locator_params(axis='x',nbins=5)
    
    # labels are identical for all plots, so just load them from the last entry in the for-loop
    handles, labels = ax.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc=(0.77, 0.20)) #TODO this might not work-well when the number of bands changes, we need to testit
    fig.supxlabel(" Magnitude ")
    fig.supylabel(" Count ")
    fig.savefig(output_directory + "mag_histogram_completeness_test.png",dpi=720,bbox_inches='tight')
    return
    
def plot_snr_mag_hist(star_catalog,gal_catalog,output_directory,SN_bins=[(9,11),(19,21)],mag_bins=np.linspace(22,26,81),skymap_center=(24000,24000),max_sep=13500):
    '''
    Draw histograms in magnitude for specific SN-slices.
    
    Args:
        star_catalog: Astropy Table; an Astropy table containing a ZP-corrected catalog of stars
        gal_catalog: Astropy Table; an Astropy table containing a ZP-corrected catalog of stars
        output_directory: string; a string specifying an output directory
        SN_bins: array-like: an array of lists specifying the range of SN-ratios to plot, defaults 9-11 (10) and 19-21 (20)
        mag_bins: Numpy array: an (N+1,) array specifying the edges of N magnitude bins
        skymap_center: list-like; a list specifying the px (x,y)-coordinates to use for a distance-cut w. max_sep
        count_range: list-like; a list specifying the y-lim of the output plot
        max_sep: float; the maximum separation for a source to be included in these plots
      
    Returns:
        None?
    
    '''
    
    col_names = star_catalog.colnames
    
    # lazy-way of getting the number of bands which exist in the catalog
    num_bands = 0
    bands_in_catalog = []
    for i in possible_bands: 
        if (i + '_psf_mag') in col_names:
            num_bands += 1
            bands_in_catalog.append(i)
    
    # cut to the center ~1deg
    #TODO center of the skymap isn't always aligned with the center of the cluster/our-pointing (although it usually is...), is there a way we can fix this even for targets where we have a non-standard pointing?
    star_center_cut = (np.sqrt((star_catalog['x'] - skymap_center[0])**2 + (star_catalog['y']-skymap_center[1])**2) < max_sep)
    gal_center_cut = (np.sqrt((gal_catalog['x'] - skymap_center[0])**2 + (gal_catalog['y']-skymap_center[1])**2) < max_sep)
    
    star_center = star_catalog[star_center_cut]
    gal_center = gal_catalog[gal_center_cut]
    
    # now create the figure and start plotting!
    rows = math.ceil((num_bands)/2)
    fig, axs = pl.subplots(2, rows, figsize=(rows*3.7,7))
    axs[-1,-1].axis('off') # clear a space for the legend later!
    
    # create a table for saving statistics
    magnitude_limit_snr = Table()
    
    # and arrays for storing the statistics
    band_array = []
    median_mag_array = []
    max_mag_array = []
    snr_array = []
    
    for i,band in enumerate(bands_in_catalog):
        
        # load the correct axes and set-up linestyle cycle
        ax = axs[i//rows, i%rows]
        styles = cycle(['-','--','-.',':'])
        
        # collect the cmodel/psf-mags, filter NaN/infs
        star_mag = star_center[band + '_psf_mag']
        gal_mag = gal_center[band + '_cmodel_mag']
        star_snr = (1/star_center[band + '_psf_magerr']) * (2.5/np.log(10))
        gal_snr = (1/gal_center[band + '_cmodel_magerr']) * (2.5/np.log(10))
        
        finite_cut_star = np.isfinite(star_mag)
        finite_cut_gal = np.isfinite(gal_mag)
        star_mag = star_mag[finite_cut_star]
        gal_mag = gal_mag[finite_cut_gal]
        star_snr = star_snr[finite_cut_star]
        gal_snr = gal_snr[finite_cut_gal]
        
        # for each SN-range, draw histogram of star/gal-SN 
        for i in range(len(SN_bins)):
            
            linestyle = next(styles)
            midbin = (SN_bins[i][0] + SN_bins[i][1])/2
            select_stars = (star_snr > SN_bins[i][0]) & (star_snr < SN_bins[i][1])
            select_gals = (gal_snr > SN_bins[i][0]) & (gal_snr < SN_bins[i][1])
            
            # stars-first, draw w. C0
            # z-band tends to be much lower-mag, so 
            if band=='z':
                mag_bin_plot = mag_bins - 1.5
            else:
                mag_bin_plot = mag_bins
            
            ax.hist(
                    star_mag[select_stars],
                    alpha = 0.9,
                    histtype='step',
                    color='C0',
                    density=True,
                    bins=mag_bin_plot,
                    label=r"Star PSF mag S/N=%s"%(int(midbin)),
                    linestyle=linestyle,
                   )
            
            # gals next, draw w. C1
            ax.hist(
                    gal_mag[select_gals],
                    alpha = 0.9,
                    histtype='step',
                    color='C1',
                    density=True,
                    bins=mag_bin_plot,
                    label=r"Galaxy CModel mag S/N=%s"%(int(midbin)),
                    linestyle=linestyle,
                   )
            
            # update the statistics
            band_array.append(band + '_psf_mag')
            median_mag_array.append(np.nanmedian(star_mag[select_stars]))
            max_mag_array.append(np.nanmax(star_mag[select_stars]))
            snr_array.append(midbin)
            
            band_array.append(band + '_cmodel_mag')
            median_mag_array.append(np.nanmedian(gal_mag[select_gals]))
            max_mag_array.append(np.nanmax(gal_mag[select_gals]))
            snr_array.append(midbin)
            
        # now make the plots pretty :)
        ax.set_title(band)
        
        # z-mag is MUCH fainter than the rest, so shift the histogram
        if band == 'z':
            ax.set_xlim((np.min(mag_bins) - 1.5,np.max(mag_bins) - 1.5))
        else:
            ax.set_xlim((np.min(mag_bins),np.max(mag_bins)))
        ax.locator_params(axis='y',nbins=3)
    
    # now update the figure labels
    handles, labels = ax.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc=(0.77, 0.20)) #TODO this might not work-well when the number of bands changes
    fig.supxlabel(" Magnitude ")
    fig.supylabel(" Frequency ")
    fig.savefig(output_directory + "snr_mag_hist.png",dpi=720,bbox_inches='tight')
    
    # update and save the table
    magnitude_limit_snr['band'] = band_array
    magnitude_limit_snr['sn_bin'] = snr_array
    magnitude_limit_snr['median'] = median_mag_array
    magnitude_limit_snr['max'] = max_mag_array
    magnitude_limit_snr.write(output_directory + "snr_mag_statistics_table.csv",format='ascii.csv',overwrite=True)
    
    return


def plot_mag_magerr(star_catalog,gal_catalog,output_directory,mag_range=(17,27),max_sep=13500,skymap_center=(24000,24000)):
    '''
    Draw a scatter-plot of the magnitude against the magnitude-error
    
    Args:
        star_catalog: Astropy Table; an Astropy table containing a ZP-corrected catalog of stars
        gal_catalog: Astropy Table; an Astropy table containing a ZP-corrected catalog of stars
        output_directory: string; a string specifying an output directory
        mag_range: array-like: an array specifying the maximum and minimum magnitude to plot
        skymap_center: list-like; a list specifying the px (x,y)-coordinates to use for a distance-cut w. max_sep
        max_sep: float; the maximum separation for a source to be included in these plots
      
    Returns:
        None?
    
    '''

    col_names = star_catalog.colnames
    
    # lazy-way of getting the number of bands which exist in the catalog
    num_bands = 0
    bands_in_catalog = []
    for i in possible_bands: 
        if (i + '_psf_mag') in col_names:
            num_bands += 1
            bands_in_catalog.append(i)
    
    # cut to the center ~1deg
    #TODO center of the skymap isn't always aligned with the center of the cluster/our-pointing (although it usually is...), is there a way we can fix this even for targets where we have a non-standard pointing?
    star_center_cut = (np.sqrt((star_catalog['x'] - skymap_center[0])**2 + (star_catalog['y']-skymap_center[1])**2) < max_sep)
    gal_center_cut = (np.sqrt((gal_catalog['x'] - skymap_center[0])**2 + (gal_catalog['y']-skymap_center[1])**2) < max_sep)
    
    star_center = star_catalog[star_center_cut]
    gal_center = gal_catalog[gal_center_cut]
    
    # now create the figure and start plotting!
    rows = math.ceil((num_bands)/2)
    fig, axs = pl.subplots(2, rows, figsize=(rows*3.7,7))
    axs[-1,-1].axis('off') # clear a space for the legend later!
    
    for i,band in enumerate(bands_in_catalog):
        
        # load the correct axes and set-up linestyle cycle
        ax = axs[i//rows, i%rows]
        styles = cycle(['-','--','-.',':'])
        
        # collect the cmodel/psf-mags, filter NaN/infs
        star_mag = star_center[band + '_psf_mag']
        gal_mag = gal_center[band + '_cmodel_mag']
        star_magerr = (star_center[band + '_psf_magerr'])
        gal_magerr = (gal_center[band + '_cmodel_magerr'])
        
        finite_cut_star = np.isfinite(star_mag)
        finite_cut_gal = np.isfinite(gal_mag)
        star_mag = star_mag[finite_cut_star]
        gal_mag = gal_mag[finite_cut_gal]
        star_magerr = star_magerr[finite_cut_star]
        gal_magerr = gal_magerr[finite_cut_gal]
        
        # draw magerr v. mag and SN-lines
        ax.scatter(star_mag, star_magerr, marker='.', alpha=0.02, label="PSF Mag (Stars)")
        ax.scatter(gal_mag, gal_magerr, marker='.', alpha=0.02, label="CModel Mag (Galaxies)")
        
        ax.hlines((2.5/np.log(10))/5, mag_range[0], mag_range[1], linestyles='solid', colors='k', alpha=0.5, label="S/N = 5")
        ax.hlines((2.5/np.log(10))/10, mag_range[0], mag_range[1], linestyles='solid', colors='k', alpha=0.5, label="S/N = 10")
        ax.hlines((2.5/np.log(10))/20, mag_range[0], mag_range[1], linestyles='solid', colors='k', alpha=0.5, label="S/N = 20")
        
        # now make the plots pretty :)
        ax.set_title(band)
        ax.set_xlim((mag_range[0],mag_range[1]))
        ax.set_ylim((0,0.25))
        ax.locator_params(axis='x',nbins=6)
    
    # now update the figure labels
    handles, labels = ax.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc=(0.77, 0.20)) #TODO this might not work-well when the number of bands changes
    for lh in leg.legend_handles:
        lh.set_alpha=0.5
    fig.supxlabel(" Magnitude ")
    fig.supylabel(" Uncertainty ")
    
    # and finally save the result
    fig.savefig(output_directory + "magerr_v_mag.png",dpi=720,bbox_inches='tight')
    return


#WARNING: column-names for gaia hard-coded since we always use gaia for photometry
def compare_astrometry(gaia_matched_catalog,catalog_tag,refcat_tag,output_directory,cluster_name,angular_bins=np.linspace(-0.4,0.4,81)):
    '''
    Draw a scatter-plot of the difference btwn LV and Gaia astrometry along with a histogram
    
    Args:
        gaia_matched_catalog: Astropy Table; a table containing matched LV-gaia objects
        catalog_tag: string; a tag appended to the end of columns from the LV catalogs
        refcat_tag: string; a tag appended to the end of columns from the Gaia catalogs
        output_directory; string; an output directory to write the figure to
        cluster_name; string; the name of the Cluster (queries coords via NED)
        angular_bins; Numpy array; an (N+1,) array of bins to use for angular displacements (arcseconds)
    
    Returns:
        None?
    
    '''
    
    # first load information from Ned about the cluster/target
    # really this is just for the declination of our pointings
    ned_result = Ned.query_object(cluster_name)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    zL = ned_result[0]['Redshift']
    
    # load the gaia and LV-coordinates; gaia/lv have the same ra/dec headers fortunately!
    lv_skycoord = SkyCoord(ra = gaia_matched_catalog['ra' + catalog_tag],dec = gaia_matched_catalog['dec' + catalog_tag],unit='deg',frame='icrs')
    gaia_skycoord = SkyCoord(ra = gaia_matched_catalog['ra' + refcat_tag],dec = gaia_matched_catalog['dec' + refcat_tag],unit='deg',frame='icrs')
    
    # compute the separations between these catalogs
    dec_sep = gaia_skycoord.dec.arcsec - lv_skycoord.dec.arcsec
    ra_sep = (gaia_skycoord.ra.arcsec - lv_skycoord.ra.arcsec)*np.cos(dec_cl * np.pi/180)
    sep = gaia_skycoord.separation(lv_skycoord)
    
    # now draw the figures
    fig, axs = pl.subplots(1, 2, figsize=(12,6))
    
    # scatter-plot of separations
    im = axs[0].scatter(ra_sep, dec_sep, marker='.', alpha=0.1)
    
    # set limits and draw a reference-grid
    axs[0].set_xlim([-0.4,0.4])
    axs[0].set_ylim([-0.4,0.4])
    axs[0].grid(color='k',alpha=0.5,linestyle=':')
    
    # figure-labels
    axs[0].set_xlabel(r"$\Delta\alpha\cos\left(\delta_{\rm cl}\right)$ ['']")
    axs[0].set_ylabel(r"$\Delta\delta$ ['']")
    
    # create and draw a 2D-kde
    #TODO there is probably a cleaner way of doing this
    delta_ra_max = np.max(ra_sep)
    delta_ra_min = np.min(ra_sep)
    delta_dec_max = np.max(dec_sep)
    delta_dec_min = np.min(dec_sep)
    
    X,Y = np.mgrid[delta_ra_min:delta_ra_max:200j, delta_dec_min:delta_dec_max:200j]
    positions = np.vstack([X.ravel(),Y.ravel()])
    values = np.vstack([ra_sep,dec_sep])
    kernel = gaussian_kde(values)
    z = kernel(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    axs[0].contour(X, Y, Z, colors='tab:red', levels=3, alpha=0.8)
    
    # now get a histogram drawn in the second pair of axes
    midpoints = (angular_bins[1:] + angular_bins[:-1])/2
    n1, b, p = axs[1].hist(ra_sep,bins=angular_bins,histtype='step',log=True,label=r"$\Delta\alpha\cos\left(\delta_{\rm cl}\right)$")
    n2, b, p = axs[1].hist(dec_sep,bins=angular_bins,histtype='step',log=True,label=r"$\Delta\delta$")
    
    # draw a reference gaussian, w/ std=0.02
    sigma = 0.02
    axs[1].plot(midpoints, norm.pdf(midpoints, 0, sigma) * (np.sqrt(2*np.pi)*sigma) * (np.max(n1) + np.max(n2)) * 0.5, alpha=0.5, linestyle=':')
    
    # make the plot pretty
    axs[1].set_ylim((8e-1, 1e4))
    axs[1].legend(loc="upper right")
    axs[1].set_xlabel("Difference (LV - Gaia) ['']")
    axs[1].set_ylabel("Count")
    
    # and save the figure
    fig.savefig(output_directory + "astrometry_check_gaia.png",dpi=720,bbox_inches='tight')
    return


def star_gal_correlation(star_catalog,gal_catalog,output_directory,cmodel_mag_cut=23.5,center_cut=13500,min_sep = 8,max_sep = 13500,skymap_center=(24000,24000),resolution=0.263,nbins=20):
    '''
    Check for correlations between the shapes of psf-stars and the measured shapes of galaxies. If our shape-measurements are correct (and remove the influence of the psf), these should be uncorrelated.
    
    Args:
        star_catalog: Astropy Table; a table of stars used for psf-measurement
        gal_catalog; Astropy Table; a table of galaxies
        output_directory: string; a filepath to an output directory
        cmodel_mag_cut: float; an upper-bound on the magnitude of galaxies used for the correlation
        center_cut: float; the maximum separation between the skymap_center and an object used for computing correlations
        min_sep: flaot; the minimum separation btwn psf-star and galaxies used for computing the correlation
        max_sep: float; the maximum separation btwn a psf-star and galaxy used for computing the correlation
        skymap_center: list; the center of the skymap
        resolution: float; resolution of the instrument ("/px)
        nbins: int; number of pins to compute the correlation in
    
    Returns:
        None?
    
    '''
    
    # cut to the center ~1deg
    #TODO center of the skymap isn't always aligned with the center of the cluster/our-pointing (although it usually is...), is there a way we can fix this even for targets where we have a non-standard pointing?
    star_center_cut = (np.sqrt((star_catalog['x'] - skymap_center[0])**2 + (star_catalog['y']-skymap_center[1])**2) < center_cut)
    gal_center_cut = (np.sqrt((gal_catalog['x'] - skymap_center[0])**2 + (gal_catalog['y']-skymap_center[1])**2) < center_cut)
    
    star_center = star_catalog[star_center_cut]
    gal_center = gal_catalog[gal_center_cut]
    
    # apply mag+lensing cut to the galaxy-catalog and select psf-stars
    gal_center_cut = gal_center[gal_center['r_cmodel_mag'] < cmodel_mag_cut]
    gal_center_cut = load_quality_cuts(gal_center_cut,quality_cuts=quality_cuts)
    star_center_cut = star_center[star_center["psf_used"] == 1]
    
    # set everythign up to compute the correlations
    cat_star = treecorr.Catalog(x=star_center_cut['x'], y=star_center_cut['y'], g1=star_center_cut['sdss_e1'], g2=star_center_cut['sdss_e2'])
    cat_galX = treecorr.Catalog(x=gal_center_cut['x'], y=gal_center_cut['y'], g1=gal_center_cut['e2'], g2=gal_center_cut['e2'])
    cat_gal10 = treecorr.Catalog(x=gal_center_cut['x'], y=gal_center_cut['y'], g1=np.zeros(len(gal_center_cut)), g2=gal_center_cut['e2'])
    cat_gal20 = treecorr.Catalog(x=gal_center_cut['x'], y=gal_center_cut['y'], g1=gal_center_cut['e1'], g2=np.zeros(len(gal_center_cut)))
    
    # compute correlations in these pixel-bins
    bin_array = np.linspace(min_sep,max_sep,nbins+1)
    midpoints = (bin_array[1:] + bin_array[:-1])/2
    
    # Finally we can compute the correlations
    # First C1, correlation btwn e2_psf and e2_gal
    gg = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_type='Linear')
    gg.process(cat_star, cat_gal20)
    C1 = gg.xip
    C1_err = gg.varxip**(0.5)

    # Next C2, correlation btwn e1_psf and e1_gal
    gg = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_type='Linear')
    gg.process(cat_star, cat_gal10)
    C2 = gg.xip
    C2_err = gg.varxip**(0.5)

    # Next C3, cross-correlation btwn e1_psf e2_gal / e2_psf e1_gal
    gg = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_type='Linear')
    gg.process(cat_star, cat_gal10)
    C3 = gg.xip
    C3_err = gg.varxip**(0.5)
    
    # In-case we need it, save a table storing the correlations
    correlation_table = Table()
    correlation_table['midpoints'] = midpoints
    correlation_table['C1'] = C1
    correlation_table['C1_err'] = C1_err
    correlation_table['C2'] = C2
    correlation_table['C2_err'] = C2_err
    correlation_table['C3'] = C3
    correlation_table['C3_err'] = C3_err
    
    correlation_table.write(output_directory + 'correlation_table.csv', format='ascii.csv', overwrite=True)
    
    # and finally draw a figure with the correlations
    pl.figure()
    
    pl.errorbar(midpoints*resolution/60, C1, yerr=C1_err, fmt='o', label='C1', alpha=0.6, capsize=3)
    pl.errorbar(midpoints*resolution/60, C2, yerr=C2_err, fmt='^', label='C2', alpha=0.6, capsize=3)
    pl.errorbar(midpoints*resolution/60, C3, yerr=C3_err, fmt='*', label='C3', alpha=0.6, capsize=3)
    pl.hlines(0, np.min(midpoints)*resolution/60, np.max(midpoints)*resolution/60, ls=':', color='k', alpha=0.5)
    
    # make the plot pretty :)
    pl.legend()
    pl.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    pl.legend()
    pl.xlabel("Distance [\']")
    pl.ylabel("Correlation")
    
    # finally save the plot
    pl.savefig(output_directory + 'correlation_plot.png', dpi=720)
    return


def draw_psf(star_catalog,output_directory,resolution=0.263,max_sep = 13500,skymap_center=(24000,24000)):
    '''
    Render shapes of the psf-stars across the coadd.
    
    Args:
        star_catalog: Astropy Table; a table of stars used for psf-measurement
        output_directory: string; a filepath to an output directory
        max_sep: float; the maximum separation btwn a psf-star and galaxy used for computing the correlation
        skymap_center: list; the center of the skymap
        resolution: float; resolution of the instrument ("/px)
    
    Returns:
        None?
    
    '''
    
    # select stars used to measure the psf
    select_psf_stars = (star_catalog['psf_used'] == 'True') & (np.sqrt((star_catalog['x'] - skymap_center[0])**2 + (star_catalog['y']-skymap_center[1])**2) < max_sep)
    psf_stars = star_catalog[select_psf_stars]
    
    col_names = star_catalog.colnames
    
    # lazy-way of getting the number of bands which exist in the catalog
    num_bands = 0
    bands_in_catalog = []
    for i in possible_bands: 
        if (i + '_psf_mag') in col_names:
            num_bands += 1
            bands_in_catalog.append(i)
    
    
    for band in ('r'): # bands_in_catalog:
        
        #kron_radius = psf_stars['rkron']
        # new export includes sdss moments, so we can compute a fwhm
        # we assume a moffat w. gamma = 2 and alpha = 1; can use these assumptions to compute the FWHM from moments
        sigma = np.sqrt(psf_stars['sdss_xx'] + psf_stars['sdss_yy']) # equivalent radius/'sigma'
        fwhm = 2*np.sqrt(np.log(2)) * sigma
        e1 = psf_stars['sdss_e1']
        e2 = psf_stars['sdss_e2']
        
        ra = psf_stars['ra']
        dec = psf_stars['dec']
        
        # plot select-visit-esque renders of the psf
        fig,ax = pl.subplots()
        im = ax.scatter(ra,dec,marker='.',c=fwhm*resolution,cmap='jet',alpha=0.4,vmin=0.7,vmax=1.5)
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("DEC (deg)")
        fig.colorbar(im)
        ax.set_title("Coadd FWHM [\"]")
        fig.savefig(output_directory + "coadd_fwhm_radius.png",dpi=720,bbox_inches='tight')
        pl.close()
        
        # draw quiver-plot showing ellipticity
        theta = np.arctan( 2 * e2/e1)/2
        ellip = np.sqrt(e1**2 + e2**2)
        u = ellip*np.cos(theta)
        v = ellip*np.sin(theta)
        
        fig,ax = pl.subplots()
        im = ax.quiver(ra,dec,u,v,ellip,cmap='jet',alpha=0.4)
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('DEC (deg)')
        fig.colorbar(im)
        ax.set_title('Coadd Ellipticity')
        fig.savefig(output_directory + "coadd_ellip.png",dpi=720,bbox_inches='tight')
        pl.close()
        
        # draw a histogram of the fwhm
        fig,ax = pl.subplots()
        median_fwhm_radius = np.nanmedian(fwhm)*resolution
        mean_fwhm_radius = np.nanmean(fwhm)*resolution
        std_fwhm_radius = np.nanstd(fwhm)*resolution
        
        ax.hist(fwhm*resolution,bins='auto',histtype='step')
        ax.set_xlabel('Seeing [\"]')
        ax.set_ylabel('Count')
        ax.set_title('Seeing Histogram')
        ax.text(0.7,0.85,r"Median FWHM: $%.2f$"%(median_fwhm_radius) + "\n" + r"Mean FWHM: $%.2f$"%(mean_fwhm_radius) + "\n" + r"Std FWHM: $%.2f$"%(std_fwhm_radius),transform=ax.transAxes)
        ax.set_xlim((0.7,1.5))
        fig.savefig(output_directory + "fwhm_hist.png",dpi=720,bbox_inches='tight')
        pl.close()
        
        # draw a histogram of the ellip
        fig,ax = pl.subplots()
        ellip = np.sqrt(e1**2 + e2**2)
        median_ellip = np.nanmedian(ellip)
        mean_ellip = np.nanmean(ellip)
        std_ellip = np.nanstd(ellip)
        
        ax.hist(ellip,bins='auto',histtype='step')
        ax.set_xlabel('Ellipticity')
        ax.set_ylabel('Count')
        ax.set_title('Ellipticity Histogram')
        ax.text(0.7,0.85,r"Median Ellip: $%.2f$"%(median_ellip) + "\n" + r"Mean Ellip: $%.2f$"%(mean_ellip) + "\n" + r"Std Ellip: $%.2f$"%(std_ellip),transform=ax.transAxes)
        fig.savefig(output_directory + "ellip_hist.png",dpi=720,bbox_inches='tight')
        pl.close()
        
        #TODO draw fwhm/moments against magnitude, good systematic check especially for BF
        
    return


def draw_visits(query_catalog,output_directory,cluster_name,instr_diameter=2.2,fov_cut=1.5):
    '''
    Draw the visits in each band.
    
    Args:
        query_catalog; Astropy Table; a table storing the queried-results from Noirlab
        output_directory: string; a filepath to an output directory
        cluster_name: string; a string specifying the name of the cluster
        instr_diameter: float; the diameter of the instrument used for imaging (deg)
        fov_cut: float; the cut applied to the fov during check/select-visit (deg)
            
    Returns:
        None?
    
    '''
    
    col_names = star_catalog.colnames
    
    # lazy-way of getting the number of bands which exist in the catalog
    num_bands = 0
    bands_in_catalog = []
    for i in possible_bands: 
        if (i + '_psf_mag') in col_names:
            num_bands += 1
            bands_in_catalog.append(i)
    
    # first load information from Ned about the cluster/target
    ned_result = Ned.query_object(cluster_name)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    zL = ned_result[0]['Redshift']
    
    # next load the coordinates of each pointing
    ra = query_catalog['ra_center']
    dec = query_catalog['dec_center']
    exposure = query_catalog['exposure']
    scale = 1/np.cos(dec_cl*np.pi/180)
    
    colors = cycle(['b','g','r','orange','brown'])
    
    
    # create figures for each band showcasing the coverage
    for i,band in enumerate(bands_in_catalog):
        
        draw_color = next(colors)
        fig,ax = pl.subplots(figsize=(6,6))
        select_band = query_catalog['ifilter'] == filter_map[band]
        ra_plot = ra[select_band]
        dec_plot = dec[select_band]
        exposure_plot = exposure[select_band]
        
        # draw a circle covering DECam's FOV at each ra/dec
        for i in range(len(ra_plot)):
            c = Ellipse((ra_plot[i],dec_plot[i]),instr_diameter*scale,instr_diameter,edgecolor='none',facecolor=draw_color,alpha=float(exposure_plot[i]/2e3))
            ax.add_patch(c)
        
        # draw a circle for our distance-cut in check_visit
        c = Ellipse((ra_cl,dec_cl),2*fov_cut*scale,2*fov_cut,edgecolor='black',facecolor='none',alpha=0.6,linestyle='--')
        ax.add_patch(c)
        
        # lastly draw a marker located at the cluster
        ax.scatter(ra_cl,dec_cl,marker='*',s=5,color='black')
        ax.set_xlabel(' RA [deg] ')
        ax.set_ylabel(' DEC [deg] ')
        
        ax.set_xlim((ra_cl + 2*scale,ra_cl - 2*scale))
        ax.set_ylim((dec_cl - 2, dec_cl + 2))
        
        # and finally save the figure
        fig.savefig(output_directory + f'visits_{band}.png',dpi=720,bbox_inches='tight')
    
    return

def draw_catalog(star_catalog,gal_catalog,output_directory,cluster_name,sn_cut=10,mag_lims=[17,27]):
    '''
    Draw the entire catalog across the FoV.
    
    Args:
        star_catalog: Astropy Table; a table of stars used for psf-measurement
        gal_catalog; Astropy Table; a table of galaxies
        output_directory: string; a filepath to an output directory
        cluster_name: string; the name of the cluster
        sn_cut: The minimum-sn of a source rendered
        mag_lims: flaot; the range of magnitudes to plot

    Returns:
        None?
    
    '''
    
    col_names = star_catalog.colnames
    
    # lazy-way of getting the number of bands which exist in the catalog
    num_bands = 0
    bands_in_catalog = []
    for i in possible_bands: 
        if (i + '_psf_mag') in col_names:
            num_bands += 1
            bands_in_catalog.append(i)
    
    # load information from Ned about the cluster/target
    ned_result = Ned.query_object(cluster_name)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    zL = ned_result[0]['Redshift']
    
    #TODO if we tweak the star/gal separation and have failed-classifications, we might drop objects we'd like to otherwise include by re-assembling the all-catalog
    all_center = vstack((star_catalog,gal_catalog))
    
    # now create the figure and start plotting!
    rows = math.ceil((num_bands)/2)
    fig, axs = pl.subplots(2, rows, figsize=(rows*3.7,7))
    axs[-1,-1].axis('off') # clear a space for the legend later!
    
    for i,band in enumerate(bands_in_catalog):
        
        # load the correct axes
        ax = axs[i//rows, i%rows]
        
        # apply the sn-cut
        select_sn = all_center[band + '_cmodel_magerr'] < (1/sn_cut)*(2.5/np.log(10))
        cmodel_mag = all_center[select_sn][band + '_cmodel_mag']
        ra = all_center[select_sn]['ra']
        dec = all_center[select_sn]['dec']
        
        im = ax.scatter(ra,dec,c=cmodel_mag,vmin=mag_lims[0],vmax=mag_lims[1],s=1,alpha=0.1)
        #cb = pl.colorbar(im,ax=ax)
        
        ax.set_title('%s-S/N > 10: %d objects'%(band, len(cmodel_mag)))
    
        # draw a marker located at the cluster
        ax.scatter(ra_cl,dec_cl,s=5,color='black',marker='*')
    
    #TODO there is a way of drawing a colorbar without changing the aspect-ratio of the axes that we should implement here
    cb = fig.colorbar(im,ax=axs.ravel().tolist())
    cb.solids.set(alpha=1)
    
    fig.supxlabel(' RA [deg] ')
    fig.supylabel(' DEC [deg] ')
    fig.savefig(output_directory + 'catalog_render.png',dpi=720,bbox_inches='tight')
    
    return

if __name__ == '__main__':
    
    # collecting arguments from cln
    if len(sys.argv) != 10:
        raise Exception("Improper usage! Correct usage: python star_catalog_filename gal_catalog_filename gaia_matched_catalog_filename catalog_tag refcat_tag query_results_filename output_directory instrument cluster_name")
    
    star_catalog_filename = sys.argv[1]
    gal_catalog_filename = sys.argv[2]
    gaia_matched_catalog_filename = sys.argv[3]
    catalog_tag = sys.argv[4]
    refcat_tag = sys.argv[5]
    query_results_filename = sys.argv[6]
    output_directory = sys.argv[7]
    instrument = sys.argv[8]
    cluster_name = sys.argv[9]
    
    # load files and thingz
    query_catalog = ascii.read(query_results_filename,format='csv')
    star_catalog = ascii.read(star_catalog_filename)
    gal_catalog = ascii.read(gal_catalog_filename)
    gaia_matched_catalog = ascii.read(gaia_matched_catalog_filename)
    resolution = instrument_resolution[instrument]
    
    # draw figures for completeness-tests
    plot_mag_hist(star_catalog,gal_catalog,output_directory)
    
    # plot histgrams in different SN-bins
    plot_snr_mag_hist(star_catalog,gal_catalog,output_directory)
    
    # draw mag v. magerr
    plot_mag_magerr(star_catalog,gal_catalog,output_directory)
    
    # check the astrometry relative to gaia
    compare_astrometry(gaia_matched_catalog,catalog_tag,refcat_tag,output_directory,cluster_name)
    
    # star-gal correlation
    star_gal_correlation(star_catalog,gal_catalog,output_directory,resolution=resolution)
    
    # draw the psf
    draw_psf(star_catalog,output_directory,resolution=resolution)
    
    # draw visits
    draw_visits(query_catalog,output_directory,cluster_name)
    
    # draw the catalog 
    draw_catalog(star_catalog,gal_catalog,output_directory,cluster_name)
    

