# actually run the photo_z step!

import sys
import os

from multiprocessing import Pool

from importlib import resources as impresources
from pathlib import Path

flux_db = impresources.files('flux_cache')

import numpy as np

from scipy.integrate import cumulative_trapezoid as cumtrapz

import matplotlib.pyplot as pl

from astropy.io import ascii
from astropy.table import Table, hstack, vstack
from astropy.cosmology import FlatLambdaCDM
from astropy.convolution import convolve, Gaussian1DKernel

# homebrew-imports here
# our default list of templates
from sed_template import default_sed_list
from filter_transmissions import DECamTransmission, SDSSTransmission, HSCTransmission

#TODO overhaul this to a standard-IO
def get_filter_map():
    '''
    This is a temporary function which outputs a hand-coded dictionary of column names mapping our catalog magnitudes to the correct header according to the filter_transmission. This requires no arguments, but in the future needs to be read from a config file
    
    Args:
        None
    
    Returns:
        filter_map: Dictionary; a dictionary which maps magnitudes to the appropriate filter in the transmission-object
        error_tag: string; tag to append to a column to get the corresponding magnitude-error

    '''
    
    
    filter_map = {
                  'u':'u_cmodel_mag',
                  'g':'g_cmodel_mag',
                  'r':'r_cmodel_mag',
                  'i':'i_cmodel_mag',
                  'z':'z_cmodel_mag',
                 }
                 
    error_tag = 'err'
    
    return filter_map, error_tag


# compute F00,F0T,FTT (Benitez+00 Eq. 9)
def compute_F(sed_template,magnitude_table,filter_map,error_tag,sed_template_upper=None,weight=None,filter_template_zps=None):
    '''
    
    Args:
        sed_template: TemplateBase; an sed_template which stores magnitudes as a function of redshifts for a given transmission-system
        magnitude_table: Astropy Table; a table containing mangitudes and magerr for all N objects you want to compute photo-z's
        filter_map: dictionary; a dictionary storing mappings from the filter-transmission to the appropriate column of magnitude_table
        error_tag: string; a tag to append to a column name to get the corresponding error
        sed_template_upper: TemplateBase; an sed_template sotirng magnitude as a function of redshift for a given transmission-system. If provided, flux from sed_template and sed_template_upper are scaled by weight and 1-weght, then combined
        weight: float; a weight to use for scaling template_upper 
        filter_template_zps: Dictionary; a dict of effective "zero-points" for the template filters
        
    Returns:
        F00: numpy array; an array of length-N for storing F00
        F0T: numpy array; an N x len(redshifts) numpy array storing F0T
        FTT: numpy array; an N x len(redshifts) numpy array storing FTT
        interpolated_scaled_flux: numpy array; an N x len(redshifts) x len(filter_map) numpy array storing the scaled template-flux
    '''
    
    # check if the sed is cached, compute it if not
    if sed_template._flux_cache == None:
        #print('Now computing fluxes for %s'%(sed_template.sed_name))
        template_flux_table = sed_template.compute_flux()
    else:
        #print("Pulling fluxes from the \'cache\'")
        template_flux_table = sed_template._flux_cache
    
    if (sed_template_upper._flux_cache == None) & (sed_template_upper is not None):
        #print('Now computing fluxes for %s'%(sed_template_upper.sed_name))
        template_flux_table_upper = sed_template_upper.compute_flux()
    elif sed_template_upper is not None:
        #print("Pulling fluxes from the \'cache\'")
        template_flux_table_upper = sed_template_upper._flux_cache
    
    # load the transmission system
    transmission = sed_template.filter_transmissions
    
    #TODO can probably do some fancy numpy tricks to compute this over all filters...
    F00 = np.zeros((len(magnitude_table),len(filter_map)))
    FTT = np.zeros((len(magnitude_table),len(sed_template.redshifts),len(filter_map)))
    F0T = np.zeros((len(magnitude_table),len(sed_template.redshifts),len(filter_map)))
    interpolated_scaled_flux = np.zeros((len(magnitude_table),len(sed_template.redshifts),len(filter_map)))
    for i in range(len(filter_map)):
        
        transmission_filter = list(filter_map.keys())[i]
        column_name = filter_map[transmission_filter]
        template_flux = template_flux_table[transmission_filter]
        if sed_template_upper is not None:
            template_flux_upper = template_flux_table_upper[transmission_filter]
            template_flux = (1-weight)*template_flux + (weight)*template_flux_upper
        
        if filter_template_zps is not None:
            template_flux = template_flux * (10**(0.4*filter_template_zps[transmission_filter]))
        
        mag_error = magnitude_table[column_name + error_tag]
        mag = magnitude_table[column_name]
        flux = 10**(-0.4 * mag)
        flux_error = ( np.log(10)/2.5 ) * np.sqrt(mag_error**2 + 0.002**2) * flux
        F00[:,i] = (flux**2) / (flux_error**2)
        F0T[:,:,i] = (template_flux * flux[:,None]) / (flux_error[:,None]**2)
        FTT[:,:,i] = (template_flux**2) / (flux_error[:,None]**2)
        interpolated_scaled_flux[:,:,i] = template_flux
        
    # this implicitly sums over all bands where the photometry exists
    F00 = np.nansum(F00,axis=1)
    F0T = np.nansum(F0T,axis=2)
    FTT = np.nansum(FTT,axis=2)
    interpolated_scaled_flux = interpolated_scaled_flux * (F0T/FTT)[:,:,None]
    
    return F00, F0T, FTT, interpolated_scaled_flux


def compute_statistics(table,template_list,trans,filter_map,error_tag,cluster_prior=None,cluster_prior_weights=None,z_spec_header=None,filter_template_zps=None,interp_steps=2,reference_mag='i',cosmo=None,redshifts=np.arange(0.01,1.5,0.001),estimator='mode',delta=0.1):
    '''
    
    Args:
        table: Astropy Table; a table containing mangitudes and magerr for all objects you want to compute photo-z's for with headers specified by filter_map
        template_list: array; an array of TemplateBase objects specifying the base templates
        trans: TransmissionBase; a class specifying the transmission system of an instrument
        filter_map: dictionary; a dictionary storing mappings from the filter-transmission to the appropriate column of magnitude_table
        cluster_prior: Numpy array; a len(redshifts) x sed_len containing an additional prior to include, weighted by cluster_prior_weight
        cluster_prior_weights: Numpy array; a 1 x sed_len array containing weights representing the fraction of galaxies in the cluster (relative to other galaxies at this same redshift)
        error_tag: string; a tag to append to a column name to get the corresponding error
        z_spec_header: string; a string specifying the column containing spec-z's. Optional, used for computing additional statistics
        filter_template_zps: Dictionary; a dict containing effective 'zero-points' for each of the filter templates, keyed by filter_map
        interp_steps: int; an integer specifying the number of interpolations to make between templates, default 2
        reference_mag: string; a string specifying the reference magnitude m0, default 'i'
        cosmo: Astropy Cosmology; an optional cosmology passed to the template_list
        redshifts: numpy array; an array storing the redshifts used to compute p(z|C,m0); defaults to np.arange(0.01,1.5,0.001)
        primary_estimate: string; a string specifying whether to use the median/mode as the primary bayesian-estimate; defaults to 'mode'
        delta: float; default delta to use for computing ODDS; defaults to 0.1
        
    Returns:
        redshift_statistics: Astropy Table; a table containing the median/mode redshift estimates, ODDS, the modified chi2, and the interpolated-index of the best-fit template
        distributions: numpy array; an an N x len(redshifts) numpy array storing p(z|C,m0)
    '''
    
    # first compute the total number of templates (including interp_steps)
    sed_len = len(template_list)
    total_templates = sed_len + (sed_len - 1)*interp_steps
    
    # initialize an array for storing the final statistics
    redshift_statistics = Table()
    
    # initialize arrays for storing intermediates
    chi2 = np.zeros((len(table),len(redshifts),(total_templates)))
    FTT_arr = np.zeros((len(table),len(redshifts),(total_templates)))
    template_prior = np.zeros((len(table),(total_templates)))
    redshift_prior = np.zeros((len(table),len(redshifts),(total_templates)))
    
    # an array used for storing the scaled flux per template
    flux_per_template = np.zeros((len(table),len(redshifts),len(filter_map),(total_templates)))
    
    for i in range(total_templates):
        
        lower_dex = int(np.floor( i / (interp_steps + 1) ))
        upper_dex = int(np.ceil( i / (interp_steps + 1) ))
        weight = ( i/ (interp_steps + 1) ) - lower_dex
        
        sed = default_sed_list[lower_dex](trans,redshifts=redshifts,cosmology=cosmo)
        sed_upper = default_sed_list[upper_dex](trans,redshifts=redshifts,cosmology=cosmo)
        name = sed.sed_name
        name_upper = sed_upper.sed_name
        
        #TODO properly search for a cache elsewhere on-disk, this will do for now
        #TODO point to flux_cache folder in the package already (somehow...)
        sed_cache = flux_db.joinpath(f'{name}_flux_cache.csv')
        sed_cache_upper = flux_db.joinpath(f'{name_upper}_flux_cache.csv')
        
        if os.path.exists(sed_cache):
            sed.load_flux_from_disk(sed_cache,trans)
            sed_upper.load_flux_from_disk(sed_cache_upper,trans)
        else:
            sed.write_flux_to_disk(sed_cache)
        
        F00,FT0,FTT,interpolated_scaled_flux = compute_F(sed,table,filter_map,error_tag,sed_upper,weight,filter_template_zps=filter_template_zps)
        chi2T = F00[:,None] - ( (FT0**2)/(FTT) )
        chi2[:,:,i] = chi2T
        FTT_arr[:,:,i] = FTT
        
        #TODO is a linear-interpolation approrpriate, or should we be doing something more sophisticated here (e.g. splines, etc.)?
        # linear interpolation for the priors as well
        template_prior[:,i] = (1-weight)*sed.get_template_prior(table[filter_map[reference_mag]]) + weight*sed_upper.get_template_prior(table[filter_map[reference_mag]])
        redshift_prior[:,:,i] = (1-weight)*sed.get_redshift_prior(table[filter_map[reference_mag]]) + weight*sed_upper.get_redshift_prior(table[filter_map[reference_mag]])
                
        # add an additional distribution to the redshift-prior to account for a cluster in the field
        if (cluster_prior != None) & (cluster_prior_weights != None):
            
            # first interpolate between templates
            interpolated_cluster_prior = (1-weight)*cluster_prior[:,lower_dex] + (weight)*cluster_prior[:,upper_dex]
            interpolated_cluster_prior_weights = (1-weight)*cluster_prior_weights[lower_dex] + (weight)*cluster_prior_weights[upper_dex]
            
            # now add this to the redshift_prior
            previous_prior = redshift_prior[:,:,i]
            redshift_prior[:,:,i] = (1 - interpolated_cluster_prior_weights)*previous_prior + (interpolated_cluster_prior_weights*interpolated_cluster_prior)[None,:]
        
        flux_per_template[:,:,:,i] = interpolated_scaled_flux
    
    # how should I go about catching nan-entries/rows and outputting nan for z_phot?
    
    # let's re-scale these according to the minimum chi-square
    chi2_zero = np.nanmin(np.nanmin(chi2,axis=2),axis=1)
    chi2 = chi2 - chi2_zero[:,None,None]
    
    # now I can assemble the probabilities
    prob_per_template = FTT_arr**(-0.5)*template_prior[:,None,:]*redshift_prior*(np.exp(-chi2/2))
    prob = np.nansum( FTT_arr**(-0.5)*template_prior[:,None,:]*redshift_prior*(np.exp(-chi2/2)),axis=2 )
    normalized_prob = prob/np.trapz(y=prob,x=redshifts,axis=1)[:,None]
    
    #TODO test convolving with a Gaussian, this was used in the older BPZ versions but my first test here I didn't see any difference
    #for i in range(len(table)):
    #    normalized_prob[i,:] = convolve(normalized_prob[i,:],Gaussian1DKernel(stddev=0.03))
    
    cumulative_prob = cumtrapz(y=normalized_prob,x=redshifts,initial=0,axis=1)
    
    # first compute the median/mode; use mode as the primary statistic in-case of multiple peaks
    z_median_indices = np.argmin(np.abs(cumulative_prob - 0.5),axis=1)
    z_median = redshifts[z_median_indices]
    z_mode_indices = np.argmax(prob,axis=1)
    z_mode = redshifts[z_mode_indices]
    
    # first find max-prob along redshift-axis (axis=1), then check of those maximized probabilities, which is largest for the best-fit template
    best_template_index = np.argmax(np.nanmax(prob_per_template,axis=1),axis=1)
    
    # selecting the primary-statistic for estimating the redshift
    if estimator == 'mode':
        primary_estimate = z_mode
        primary_estimate_indices = z_mode_indices
    elif estimator == 'median':
        primary_estimate = z_median
        primary_estimate_indices = z_median_indices
    else:
        print(f'{estimator} is not one of median or mode... defaulting to mode!')
        primary_estimate = z_mode
        primary_estimate_indices = z_mode_indices
    
    # computing ODDS
    upper = np.argmin(np.abs(redshifts - (primary_estimate[:,None] - delta*(1+primary_estimate[:,None]))),axis=1)
    lower = np.argmin(np.abs(redshifts - (primary_estimate[:,None] + delta*(1+primary_estimate[:,None]))),axis=1)
    
    best_fit_flux = np.zeros((len(table),len(filter_map)))
    odds = np.zeros(len(table))
    
    #TODO for-loop here isn't great, there should be a way of doing this using broadcasting
    for i in range(len(table)):
        best_fit_flux[i] = flux_per_template[i,primary_estimate_indices[i],:,best_template_index[i]]
        odds[i] = np.abs(cumulative_prob[i,upper[i]] - cumulative_prob[i,lower[i]])
        
    template_error = np.nanmax(best_fit_flux,axis=1)/15
    
    # now go through each band to compute chi2 mod
    mod_chi2_per_band = np.zeros((len(table),len(filter_map)))
    for i in range(len(filter_map)):
        
        transmission_filter = list(filter_map.keys())[i]
        column_name = filter_map[transmission_filter]
        
        observed_flux = 10**(-0.4 * table[column_name])
        observed_flux_err = ( np.log(10)/2.5 ) * np.sqrt(table[column_name + error_tag]**2 + 0.002**2) * observed_flux
        
        mod_chi2_per_band[:,i] = (( observed_flux - best_fit_flux[:,i] )**2)/( observed_flux_err**2 + template_error**2 )
    
    # compute the number of finite fluxes to get the DOF, assume DOF=1 for bands<4
    mod_chi2_dof = np.sum(np.isfinite(mod_chi2_per_band),axis=1) - 3
    mod_chi2_dof[mod_chi2_dof < 1] = 1
    mod_chi2 = np.nansum(mod_chi2_per_band,axis=1)/mod_chi2_dof
            
    # and finally, I can wrap all the statistics into one output-table
    #TODO add ML (min-chi2) estimate to outputs
    redshift_statistics['z_phot'] = primary_estimate
    redshift_statistics['z_median'] = z_median
    redshift_statistics['z_mode'] = z_mode
    redshift_statistics['template_at_z_phot'] = best_template_index
    redshift_statistics['odds'] = odds
    redshift_statistics['mod_chi2'] = mod_chi2
    
    for i in range(len(filter_map)):
        
        transmission_filter = list(filter_map.keys())[i]
        column_name = filter_map[transmission_filter]

        redshift_statistics[f'{transmission_filter}_flux_at_z_phot'] = best_fit_flux[:,i]
    
    # if z_spec is given, run quality checks based on it
    if z_spec_header is not None:
        
        z_spec = table[z_spec_header]
        z_spec_index = np.argmin(np.abs(z_spec[:,None] - redshifts),axis=1)
        
        # create some arrays for storing outputs        
        ML_template_at_z_spec = np.zeros(len(table))
        flux_at_z_spec = np.zeros((len(table),len(filter_map)))
        
        #TODO another for-loop here, I can't think of a clever way of doing this right now
        for i in range(len(table)):
        
            template_z_spec = np.argmin(chi2[i,z_spec_index[i],:])
            
            ML_template_at_z_spec[i] = template_z_spec
            flux_at_z_spec[i,:] = flux_per_template[i,z_spec_index[i],:,template_z_spec]
        
        redshift_statistics['template_at_z_spec'] = template_z_spec
        
        for i in range(len(filter_map)):
        
            transmission_filter = list(filter_map.keys())[i]
            column_name = filter_map[transmission_filter]
            
            redshift_statistics[f'{transmission_filter}_flux_at_z_spec'] = flux_at_z_spec[:,i]

    return redshift_statistics, normalized_prob

# a helper function for binning arrays according to bins defined in some ind. variable
def bin_array_in_xbins(array,x,xbin_edges):
    '''
    
    Args:
        array: N x 1 Numpy array; an N-dim array containing the quantity to bin
        x: N x 1 Numpy array; an N-dim array specifying the value of an ind. variable for each value in array
        xbin_edges: Numpy array; an array of length M+1 specifying the edges of M-bins in x 
    Returns:
        mean: M x 1 Numpy array; the mean of array in each xbin_edges
        median: M x 1 Numpy array; the median of array in each xbin_edges
        std: M x 1 Numpy array; the stdev of array in each xbin_edges
        binned_array: N x M Numpy array; an N x M boolean array specifying whether entry (i,j) belongs in bin-j
    '''

    mean = np.zeros(len(xbin_edges - 1))
    median = np.zeros(len(xbin_edges - 1))
    std = np.zeros(len(xbin_edges - 1))
    binned_array = np.zeros((len(array),len(xbin_edges)-1),dtype=bool)
    for i in range(len(xbin_edges) - 1):
        bin_bool = (x >= xbin_edges[i]) & (x < xbin_edges[i+1])
        binned_array[:,i] = bin_bool
        mean[i] = np.nanmean(array[bin_bool])
        median[i] = np.nanmedian(array[bin_bool])
        std[i] = np.nanstd(array[bin_bool])
    
    return mean, median, std, binned_array

'''
Use this to draw plots to assess the quality of photo-zs, including the following:

1. Histogram of (zs - zp)/(1+zs) (with corresponding gaussian and marks where the outliar-cut is located)
2. zp v zs w/ lines at the boundary of catastrophic-outliars (the dashed-line things), colored by chi2_mod
3. plot chi2_mod/odds against blendedness, colored by magnitude
4. (zs - zp)/(1+zs) v blendedness, colored by magnitude
4. maybe ELP if I can find out what it means?

I'll write the code for these, but they'll only work really well with a large spec-z redshift sample (e.g. collect multiple gen3 clusters together and run bpz to see the results)

1. median-bias vs mag/redshift
2. cat-outlier vs mag/redshift 
3. stdv (outliers-rejected) vs mag/redshift

'''

def large_sample_quality_check_plots(redshift_statistics,table,filter_map,z_spec_header=None,output_filepath='',blendedness_header='blendedness',refmag='i'):
    '''
    
    Args:
        redshift_statistics: Astropy Table; a table containing the median/mode redshift estimates, ODDS, the modified chi2, and the interpolated-index of the best-fit template
        spec_z_header: string; a string conatining the header of spec_z's in table, defaults to None
        output_filepath: string; a string specifying the filepath to be appended to the beginning of the output figures
        table: Astropy Table; ; a table containing mangitudes and magerr for all objects you want to compute photo-z's for with headers specified by filter_map.
        blendedness_header: string; string specifying the column-name for the blendedness of an object
        
    Returns:
        None?

    '''
    
    # collect the arrays I'll be needing from table/redshift_statistics
    odds = redshift_statistics['odds']
    mod_chi2 = redshift_statistics['mod_chi2']
    z_phot = redshift_statistics['z_phot']
    ref_mag = table[filter_map[refmag]]
    blendedness = table[blendedness_header]
    z_spec = table[z_spec_header]
    
    # there are a few easy-to-compute statisticds I need to check
    num_photoz = np.sum(np.isfinite(z_phot))
    delta_z = (z_spec - z_phot)/(1+z_spec)
    cat_outliers = np.nansum( np.abs(delta_z) > 0.15 )/num_photoz
    cat_outliers_percent = cat_outliers*100
    NMAD = np.nanmedian( 1.48 * np.abs(delta_z) )
    

    # first, generic zs vs zp plot + colorbar for reduced chi-2
    fig,ax = pl.subplots()
    im = ax.scatter(x=z_spec,y=z_phot,s=5,c=mod_chi2,alpha=0.7,vmin=0,vmax=1,cmap='Spectral')
    cb = fig.colorbar(im)
    
    # add-in curves for 5% and 15% difference in redshift
    x_val = np.linspace(-0.1, 1.1, 10)
    ax.plot(x_val, x_val, 'k-', alpha=0.5)
    ax.plot(x_val, x_val + 0.05 * (1. + x_val), 'k--', alpha=0.5)
    ax.plot(x_val, x_val - 0.05 * (1. + x_val), 'k--', alpha=0.5)
    ax.plot(x_val, x_val + 0.15 * (1. + x_val), 'k:', alpha=0.5)
    ax.plot(x_val, x_val - 0.15 * (1. + x_val), 'k:', alpha=0.5)
    
    ax.set_xlim((-0.001,1.001))
    ax.set_ylim((-0.001,1.001))
    
    # setting labels
    cb.ax.set_ylabel('$ \\chi_{m}^2 $',rotation='horizontal',verticalalignment='center',size=12)
    ax.set_xlabel('$ z_s $')
    ax.set_ylabel('$ z_p $')
    ax.set_title(f'{num_photoz} Galaxies: NMAD={NMAD:.3f} and $ \\eta $={cat_outliers_percent:.1f}%')
    fig.savefig(output_filepath + 'photz_v_specz_chi2_color_before_qc.png',bbox_inches='tight',dpi=720)
    
    
    # plot cat-outlier-rate against magnitude/redshift
    fig,(ax1,ax2) = pl.subplots(1,2,sharey=True,gridspec_kw={'wspace':0})
    
    # plot stdv of delta_z in mag/redshift bins
    fig_2,(ax1_2,ax2_2) = pl.subplots(1,2,sharey=True,gridspec_kw={'wspace':0})
    
    # define redshift and magnitude bins
    redshift_bins = np.arange(0.01,1,0.075)
    mag_bins = np.arange(13,26,1)
    
    # now bin modified-chi2 as a function of z_spec
    mean,median,std,bin_bool = bin_array_in_xbins(delta_z,z_spec,redshift_bins)
    mean_z,median_z,std_z,bin_bool_z = bin_array_in_xbins(z_spec,z_spec,redshift_bins)
    
    # draw the std-per-bin
    non_zero_filter = ( np.isfinite(std) & (std > 0) )
    ax1_2.plot(median_z[non_zero_filter],std[non_zero_filter])
    
    # first for redshift bins
    outlier_rate_per_bin = []
    errorbars = []
    for i in range(len(redshift_bins) - 1):
        dz_per_bin = delta_z[bin_bool[:,i]]
        errorbars.append(1/np.sqrt(np.sum(np.isfinite(dz_per_bin))))
        outlier_rate_per_bin.append(np.sum(np.abs(dz_per_bin[np.isfinite(dz_per_bin)]) > 0.15)/np.sum(np.isfinite(dz_per_bin)))
    
    outlier_rate_per_zbin = np.array(outlier_rate_per_bin)
    errorbars_zbin = np.array(errorbars)
    
    # also bin modified-chi2 as a function of reg_mag
    mean,median,std,bin_bool = bin_array_in_xbins(delta_z,ref_mag,mag_bins)
    mean_mag,median_mag,std_mag,bin_bool_mag = bin_array_in_xbins(ref_mag,ref_mag,mag_bins)
    
    # draw the std-per-bin
    non_zero_filter = ( np.isfinite(std) & (std > 0) )
    ax2_2.plot(median_mag[non_zero_filter],std[non_zero_filter])
    
    # next for magnitude bins
    outlier_rate_per_bin = []
    errorbars = []
    for i in range(len(mag_bins) - 1):
        dz_per_bin = delta_z[bin_bool[:,i]]
        errorbars.append(1/np.sqrt(np.sum(np.isfinite(dz_per_bin))))
        outlier_rate_per_bin.append(np.sum(np.abs(dz_per_bin[np.isfinite(dz_per_bin)]) > 0.15)/np.sum(np.isfinite(dz_per_bin)))

    outlier_rate_per_mbin = np.array(outlier_rate_per_bin)
    errorbars_mbin = np.array(errorbars)
    
    # plotting cat-outlier rate
    ax1.fill_between(median_z[:-1],outlier_rate_per_zbin+errorbars_zbin,outlier_rate_per_zbin-errorbars_zbin,alpha=0.5)
    ax2.fill_between(median_mag[:-1],outlier_rate_per_mbin+errorbars_mbin,outlier_rate_per_mbin-errorbars_mbin,alpha=0.5) 
    
    ax1.set_xlabel(' z$_{spec}$ ')
    ax2.set_xlabel(' i ')
    ax1.set_ylabel('$ \\eta $')
    ax1.set_ylim((-0.01,0.61))
    
    fig.savefig(output_filepath + 'outlier_rate_v_mag_z_no_qc.png',bbox_inches='tight',dpi=720)
    
    # plotting std of delta_z
    ax1_2.set_xlabel(' z$_{spec}$ ')
    ax2_2.set_xlabel(' i ')
    ax1_2.set_ylabel('$ \\sigma_{\\Delta z} $')
    fig_2.savefig(output_filepath + 'delta_z_std_v_mag_redshift_no_qc.png',bbox_inches='tight',dpi=720)
    
    # next, plot eta and sigma
    fig,ax = pl.subplots()
    chi2_bins = np.arange(0,8,1)
    mean,median,mode,bin_bool = bin_array_in_xbins(delta_z,mod_chi2,chi2_bins)
    mean_chi2,median_chi2,mode_chi2,bin_bool_chi2 = bin_array_in_xbins(mod_chi2,mod_chi2,chi2_bins)
    
    # next for magnitude bins
    outlier_rate_per_bin = []
    errorbars = []
    for i in range(len(chi2_bins) - 1):
        dz_per_bin = delta_z[bin_bool[:,i]]
        errorbars.append(1/np.sqrt(np.sum(np.isfinite(dz_per_bin))))
        outlier_rate_per_bin.append(np.sum(np.abs(dz_per_bin[np.isfinite(dz_per_bin)]) > 0.15)/np.sum(np.isfinite(dz_per_bin)))

    outlier_rate_per_chi2= np.array(outlier_rate_per_bin)
    errorbars_chi2 = np.array(errorbars)
    
    ax.fill_between(median_chi2[:-1],outlier_rate_per_chi2+errorbars_chi2,outlier_rate_per_chi2-errorbars_chi2,alpha=0.5)
    ax.set_xlabel('$ \\chi_{m}^2 $')
    ax.set_ylabel('$ \\eta $')
    
    fig.savefig(output_filepath + 'outlier_rate_v_mod_chi2.png',bbox_inches='tight',dpi=720)
    
    
    #TODO draw u-g v r-i with templates for outliers, colored by z_phot, z_spec, delta_z, ref_mag
    # That will be very useful for experimenting with the addition of 'physical priors', e.g. intrinsic-extinction
    
    return


def generic_quality_check_plots(redshift_statistics,table,filter_map,z_spec_header=None,output_filepath='',blendedness_header='blendedness',refmag='i',mod_chi2_cut=4,odds_cut=0.95):
    '''
    
    Args:
        redshift_statistics: Astropy Table; a table containing the median/mode redshift estimates, ODDS, the modified chi2, and the interpolated-index of the best-fit template
        spec_z_header: string; a string conatining the header of spec_z's in table, defaults to None
        output_filepath: string; a string specifying the filepath to be appended to the beginning of the output figures
        table: Astropy Table; a table containing mangitudes and magerr for all objects you want to compute photo-z's for with headers specified by filter_map.
        blendedness_header: string; string specifying the column-name for the blendedness of an object
        
    Returns:
        None?

    '''
    
    # collect the arrays I'll be needing from table/redshift_statistics
    odds = redshift_statistics['odds']
    mod_chi2 = redshift_statistics['mod_chi2']
    
    # apply qc if given
    if mod_chi2_cut is not None:
        cut = mod_chi2 < mod_chi2_cut
        table = table[cut]
        redshift_statistics = redshift_statistics[cut]
        odds = odds[cut]
        mod_chi2 = mod_chi2[cut]
    
    if odds_cut is not None:
        cut = odds > odds_cut
        table = table[cut]
        redshift_statistics = redshift_statistics[cut]
        odds = odds[cut]
        mod_chi2 = mod_chi2[cut]
    
    # computing the remaining statistics
    z_phot = redshift_statistics['z_phot']
    ref_mag = table[filter_map[refmag]]
    blendedness = table[blendedness_header]
    z_spec = table[z_spec_header]
    
    # there are a few easy-to-compute statisticds I need to check
    num_photoz = np.sum(np.isfinite(z_phot))
    delta_z = (z_spec - z_phot)/(1+z_spec)
    cat_outliers = np.nansum( np.abs(delta_z) > 0.15 )/num_photoz
    cat_outliers_percent = cat_outliers*100
    NMAD = np.nanmedian( 1.48 * np.abs(delta_z) )
    
    
    # first, generic zs vs zp plot + colorbar for reduced chi-2
    fig,ax = pl.subplots()
    im = ax.scatter(x=z_spec,y=z_phot,s=5,c=mod_chi2,alpha=0.7,vmin=0,vmax=1,cmap='Spectral')
    cb = fig.colorbar(im)
    
    # add-in curves for 5% and 15% difference in redshift
    x_val = np.linspace(-0.1, 1.1, 10)
    ax.plot(x_val, x_val, 'k-', alpha=0.5)
    ax.plot(x_val, x_val + 0.05 * (1. + x_val), 'k--', alpha=0.5)
    ax.plot(x_val, x_val - 0.05 * (1. + x_val), 'k--', alpha=0.5)
    ax.plot(x_val, x_val + 0.15 * (1. + x_val), 'k:', alpha=0.5)
    ax.plot(x_val, x_val - 0.15 * (1. + x_val), 'k:', alpha=0.5)
    
    ax.set_xlim((-0.001,1.001))
    ax.set_ylim((-0.001,1.001))
    
    # setting labels
    cb.ax.set_ylabel('$ \\chi_{m}^2 $',rotation='horizontal',verticalalignment='center',size=12)
    ax.set_xlabel('$ z_s $')
    ax.set_ylabel('$ z_p $')
    ax.set_title(f'{num_photoz} Galaxies: NMAD={NMAD:.3f} and $ \\eta $={cat_outliers_percent:.1f}%')
    fig.savefig(output_filepath + 'photz_v_specz_chi2_color.png',bbox_inches='tight',dpi=720)

    # second, generic zs vs zp plot + colorbar for mag
    fig,ax = pl.subplots()
    im = ax.scatter(x=z_spec,y=z_phot,s=5,c=ref_mag,alpha=0.7,vmin=13,vmax=22)
    cb = fig.colorbar(im)
    
    # add-in curves for 5% and 15% difference in redshift
    x_val = np.linspace(-0.1, 1.1, 10)
    ax.plot(x_val, x_val, 'k-', alpha=0.5)
    ax.plot(x_val, x_val + 0.05 * (1. + x_val), 'k--', alpha=0.5)
    ax.plot(x_val, x_val - 0.05 * (1. + x_val), 'k--', alpha=0.5)
    ax.plot(x_val, x_val + 0.15 * (1. + x_val), 'k:', alpha=0.5)
    ax.plot(x_val, x_val - 0.15 * (1. + x_val), 'k:', alpha=0.5)
    
    ax.set_xlim((-0.001,1.001))
    ax.set_ylim((-0.001,1.001))
    
    # setting labels
    cb.ax.set_ylabel('$ i $',rotation='horizontal',verticalalignment='center',size=12)
    ax.set_xlabel('$ z_s $')
    ax.set_ylabel('$ z_p $')
    ax.set_title(f'{num_photoz} Galaxies: NMAD={NMAD:.3f} and $ \\eta $={cat_outliers_percent:.1f}%')
    fig.savefig(output_filepath + 'photz_v_specz_mag_color.png',bbox_inches='tight',dpi=720)
    
    
    # next, histograms of the spec-z and phot-z
    bin_edges = np.arange(0.01,1,0.05)
    fig,ax = pl.subplots()
    
    # draw the histograms
    ax.hist(z_phot,bins=bin_edges,histtype='step',alpha=0.75,label='Photo-z')
    ax.hist(z_spec,bins=bin_edges,histtype='step',alpha=0.75,label='Spec-z')
    ax.legend()
    ax.set_xlabel('Redshift')
    ax.set_ylabel('Count')
    ax.set_title('Redshift Distributions')
    fig.savefig(output_filepath + 'redshift_histograms.png',bbox_inches='tight',dpi=720)
    
    
    # and FINALLY, a histogram of delta-z
    fig,ax = pl.subplots()
    z_bins = np.arange(-0.5,0.5,0.05)
    
    # draw the bins
    ax.hist(delta_z,bins=100,density=True,histtype='step')
    
    # what the distribution should be
    x = np.linspace(-0.5,0.5,1000)
    gauss = (1/(NMAD*np.sqrt(2*np.pi))) * np.exp( (- x**2)/(2*(NMAD**2)))
    pl.plot(x,gauss,'--',alpha=0.78,label=r"$\mathcal{N}$~(0, NMAD)")
    
    # labelling everything
    ax.set_xlim((-0.55,0.55))
    ax.set_xlabel('$ \\Delta z $')
    ax.set_ylabel(' Frequency ')
    ax.legend()
    fig.savefig(output_filepath + 'delta_z_histogram.png',bbox_inches='tight',dpi=720)
    
    
if __name__ == '__main__':
    
    if len(sys.argv)==6:
        
        # collecting arguments from cln
        table_filename = sys.argv[1]
        output_catalog_filename = sys.argv[2]
        output_directory = sys.argv[3]
        cores = int(sys.argv[4])
        instrument = sys.argv[5]
        z_spec_header = None # if this isn't specified, skip generating output figures
        
    elif len(sys.argv)==7:
    
        # collecting arguments from cln
        table_filename = sys.argv[1]
        output_catalog_filename = sys.argv[2]
        output_directory = sys.argv[3]
        cores = int(sys.argv[4])
        instrument = sys.argv[5]
        z_spec_header = sys.argv[6]
    
    else:
        print("python bayesian_photo_z.py catalog output_catalog_filename output_directory cores instrument [OPTIONAL: z_spec_header]")
        raise Exception("Improper Usage! Correct usage: python bayesian_photo_z.py catalog output_catalog_filename output_directory cores instrument [OPTIONAL: z_spec_header]")
    
    #TODO should this dictionary come from our standard IO or be defined in filter_transmissions.py? I'm not sure, leaving it here for now
    # load the appropriate transmission system
    filter_transmissions = {'decam':DECamTransmission, 'sdss':SDSSTransmission, 'hsc':HSCTransmission}
    if instrument not in filter_transmissions:
        print(f'{instrument} is not in the dictionary of transmission systems!')
        raise Exception(f'{instrument} is not in the dictionary of transmission systems!')
    
    #TODO move these configurables to a standard IO
    mod_chi2_cut = 4
    odds_cut = 0.95
    
    trans = filter_transmissions[instrument]()
    filter_map,error_tag = get_filter_map()
    
    # loading the table
    table = ascii.read(table_filename)
    
    # split the table into chunks for running in parallel; approximately 60 per-task
    #chunk_size = np.ceil(len(table)/cores)
    #split_table = [table[int(i*chunk_size):int((i+1)*chunk_size)] for i in range(cores)]
    #print([(int(i*chunk_size),int((i+1)*chunk_size)) for i in range(cores)])
    
    # breakup the table into chunks of 60-entries
    # for 20 massive groups, processes are killed by OS
    chunk_size = 60
    num_chunks = int(np.ceil(len(table)/chunk_size))
    split_table = [table[int(i*chunk_size):int((i+1)*chunk_size)] for i in range(num_chunks)]
    
    # lazy trick for passing all the necessary arguments, define a function to wrap this
    def compute_statistics_parallel(x):
        out = compute_statistics(x,default_sed_list,trans,filter_map,error_tag,z_spec_header,interp_steps=2)
        print("Exiting pool...")
        return out
    
    with Pool(cores) as p:
        outputs = p.map(compute_statistics_parallel,split_table)
    
    for i in range(len(outputs)):
        if i == 0:
            redshift_statistics = outputs[i][0]
        else:
            redshift_statistics = vstack([redshift_statistics,outputs[i][0]])
    
    # if z_spec_header is defined, draw plots for the generic quality-checks
    if z_spec_header is not None:
        generic_quality_check_plots(redshift_statistics,table,filter_map,z_spec_header=z_spec_header,output_filepath=output_directory,mod_chi2_cut=mod_chi2_cut,odds_cut=odds_cut)
    
    # append redshift_statistics to the table and write it to disk
    updated_table = hstack([table,redshift_statistics])
    updated_table.write(output_directory + output_catalog_filename, format="ascii.csv", overwrite=True)
    
    #TODO recalibrate priors eventually, we're waiting to use DESI for low-z; but will likely need COSMOS or other fields for mag > ~21 
    
    #TODO now that we have this functionality, someone should experiment with an "sed-zp correction". With these initial values (from a large sample of spec-z across four clusters), our cat-error goes UP by ~1.5%, but with some tuning it could help? LePhare (another photo-z code) among others include corrections like this...
    #filter_template_zps = {'u':0.242,'g':0.011,'r':-0.006,'i':0.004,'z':-0.015}
    
    #TODO functionality for cluster priors has been added, we should experiment with this to see if there is any gain,
    # might need to pass coordinates to really gain from this, ~1-10% of galaxies across the full-fov are cluster-members, but all galaxies w/in R500 will be member and have a different dist of spectral types
    # Once we have the full LV sample processed + DESI, someone should try to recalibrate priors p(T|z) -> p(T|z,R) and p(z|T,m0) -> p(z|T,m0,R) where R is some distance from the cluster-center
    # We also need to decide if this requires any additional priors that should be placed on R... project for someone else!
    
        
        
        




