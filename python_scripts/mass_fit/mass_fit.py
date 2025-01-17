# run mass_fit

import sys
import os
from pathlib import Path

from multiprocessing import Pool

import numpy as np

from scipy.integrate import trapz, cumtrapz
from scipy.optimize import minimize

from scipy import interpolate
from scipy import stats


import matplotlib.pyplot as pl

from astropy.io import ascii, fits
from astropy.table import Table, hstack, vstack
from astropy.wcs import WCS

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

from photutils.segmentation import detect_sources, SourceCatalog

from astroquery.ipac.ned import Ned

# HOMEBREW MODULES BELOW
from shear_profiles import implemented_profiles, shear_nfw, sigma_crit


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
    #TODO is there a better way of formatting this information?
    #TODO should the column-header to r_cmodel_magerr be generalized (e.g. w/ standard-IO) or left to the user to specify correctly?
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
    
    return cut_table


# helper for fitting multiple-profiles
def multipeak_residuals(M200_array,zS,x_pos,y_pos,g1_observed,g2_observed,rad_per_px,zL=0.05,num_peaks=1,cluster_centers=[(0,0)],shear_profile=shear_nfw,model_parameters=[{}],cosmo=cosmo):
    '''
    A function for helping multipeak fitting, computes the chi2 given N-peaks
    
    Args:
        M200: numpy array; (num_peaks,)-shape numpy array specifying the mass of a given peak
        zL: float; float specifying the redshift of the peaks
        zS: numpy array; (N,)-shape array specifying the photo-z of objects to fit
        x_pos: numpy array; (N,)-shape numpy array specifying the x-position of objects to fit (pixels)
        y_pos: numpy array; (N,)-shape numpy array specifying the y-position of objects to fit (pixels)
        g1_observed: numpy array; (N,)-shape numpy array specifying the reduced-shears_1 to fit
        g2_observed: numpy array; (N,)-shape numpy array specifying the reduced-shears_2 to fit
        num_peaks: int; integer specifying the number of peaks used in the fitting
        cluster_centers: array-like; (num_peaks,)-shape array storing the pixel-coordinates of a cluster-center (pixels)
        shear_profile: function; a function specifying the shear profile, expects to return kappa, gamma1, and gamma2; defaults to shear_nfw
        model_parameters: array-like; (num_peaks,)-shape array containing kwargs to unpack for each peak
        cosmo: Astropy cosmology; a cosmology to pass to the shear_profile
    
    Returns:
        g1_residuals: numpy array; (N,)-shape numpy array storing the difference between the N-peak model and observed shears (g_1)
        g2_residuals: numpy array; (N,)-shape numpy array storing the difference between the N-peak model and observed shears (g_2)
        
    '''
    
    # initialize arrays to store the residuals
    gamma_1 = np.zeros(len(zS))
    gamma_2 = np.zeros(len(zS))
    kappa = np.zeros(len(zS))
    
    # shears/convgs are linear for objects at the same redshift, so loop over all peaks then compute the reduced-shear
    for i in range(num_peaks):
        
        delta_x = x_pos - cluster_centers[i][0]
        delta_y = y_pos - cluster_centers[i][1]
        kap,gam1,gam2 = shear_profile(M200_array[i],zS=zS,delta_x=delta_x*rad_per_px,delta_y=delta_y*rad_per_px,zL=zL,cosmo=cosmo,**model_parameters[i])

        gamma_1 += gam1
        gamma_2 += gam2
        kappa += kap
    
    # the theoretical reduced-shears for each object
    g1_theory = gamma_1/(1-kappa)
    g2_theory = gamma_2/(1-kappa)
    
    # now the residuals can be computed
    g1_residuals = g1_theory - g1_observed 
    g2_residuals = g2_theory - g2_observed
    
    return g1_residuals, g2_residuals


# helper function for drawing plots of the statistics
def draw_plots(masses,names,reduced_chi2,complete_sample,output_directory,resolution,cosmo=cosmo):
    '''
    A function for helping multipeak fitting, computes the chi2 given N-peaks.
    
    Args:
        masses: Numpy Array; an (Nboot,Npeak) array containing the best-fit mass estimates
        names: array; Names of the objects/peaks, passed to figure titles
        reduced_chi2: Numpy array; an (Nboot,)-array storing the reduced-chi2 statistic for each best-fit mass
        complete_sample: dict; a dictionary of kwargs passed to multipeak_residual for computing residuals at the median mass
        output_directory: String; a string specifying the output directory to save plots in
        resolution: float; the resolution of the instrument in "/px
    
    Returns:
        None?
    
    '''
    
    # arrays for storing useful/important statistics
    median_masses = []
    upper_masses = []
    lower_masses = []
    
    # for now, draw a separate plot for each peak
    #TODO we should automatically generate a corner-plot when there is more than one peak, apparently this isn't hard using the corner package?
    for i in range(1): #len(masses[0,:])):
        
        peak_mass = masses/(1e14) #[0,i]
        fig,ax = pl.subplots()
        kde = stats.gaussian_kde(peak_mass)
        mass_range = np.max(peak_mass) - np.min(peak_mass)
        
        # build an interpolated kde
        plot_mass_range = np.linspace(np.min(peak_mass) - mass_range/2,np.max(peak_mass) + mass_range/2, 2000)
        
        #TODO should we split these statistics into a separate function for readability?
        #TODO are these extra normalizations/interpolations necessary?
        pdf = kde(plot_mass_range)
        pdf = interpolate.interp1d(plot_mass_range, pdf / trapz(pdf,plot_mass_range),bounds_error=False,fill_value=0)
        inverse_cdf = interpolate.interp1d(cumtrapz(pdf(plot_mass_range), plot_mass_range, initial=0),plot_mass_range,bounds_error=False)
        
        median_mass = inverse_cdf(0.5)
        median_prob = pdf(median_mass)
        sig_upper = inverse_cdf(0.5 + 0.341)
        sig_lower = inverse_cdf(0.5 - 0.341)
        
        median_masses.append(median_mass)
        upper_masses.append(sig_upper)
        lower_masses.append(sig_lower)
        
        mean_mass = np.trapz(peak_mass * pdf(peak_mass), peak_mass)
        mean_prob = pdf(mean_mass)
        
        ax.hist(peak_mass, bins=20, density=True, histtype='step')
        ax.plot(plot_mass_range,pdf(plot_mass_range))
        
        ax.vlines(median_mass, 0, median_prob, linestyles='dashed', label = 'Median', color='k')
        #ax.vlines(mean_mass, 0, mean_prob, linestyles='dotted', label='Mean', color='k')
        
        ax.set_xlabel(r"$M_{200}$ [$10^{14}M_\odot$]")
        ax.set_ylabel("Frequency")
        peak_name = names[i]
        title_string = f'{peak_name}: ' + r"$M_{200} = %.2f_{%.2f}^{+ %.2f} \times 10^{%s}M_\odot$"%( 
                        median_mass, 
                        (sig_lower - median_mass), 
                        (sig_upper - median_mass),
                        14,
                    )
        
        # lazy trick, draw an invisible line to get some text in the legend
        chi2_tag = r"$\overline{\chi}^2_{\nu} = %.2f$"%(np.mean(reduced_chi2))
        ax.plot(np.NaN,np.NaN,'-',color='none',label=chi2_tag)
        
        ax.legend()
        ax.set_title(title_string)
        fig.savefig(output_directory + 'mass_histogram.png',dpi=720,bbox_inches='tight')
    
    # table for writing statistics to disk
    mass_statistics = Table()
    # mass_statistics['peak_index'] = ... #TODO for multiple peaks we need an index or flag to track the mass associated with each one
    mass_statistics['median_mass'] = np.array(median_masses)*1e14
    mass_statistics['mass_upper_ci'] = np.array(upper_masses)*1e14
    mass_statistics['mass_lower_ci'] = np.array(lower_masses)*1e14
    mass_statistics['reduced_chi_square'] = np.ones(len(median_masses))*np.mean(reduced_chi2)
    mass_statistics.write(output_directory + 'mass_peak_statistics.csv',format='ascii.csv',overwrite=True)
    
    # also write the individual bootstrap realizations to disk (will be useful for later)
    mass_realizations = Table()
    
    #TODO when we use multiple-peaks, we'll want to explore the parameter space more efficiently (MCMC?)
    # anyways as with the other loops here, I've hard-coded this to one for now...
    for i in range(1):
        mass_realizations['masses'] = masses
    mass_realizations.write(output_directory + 'mass_realizations.csv',format='ascii.csv',overwrite=True)
    
    
    # now draw a histogram of the reduced-chi2
    fig,ax = pl.subplots()
    
    ax.hist(reduced_chi2,bins=20,histtype='step',density=True)
    ax.set_xlabel(r"$\chi^2_{\nu}$")
    ax.set_ylabel("Frequency")
    fig.savefig(output_directory + 'chi2_histogram.png',dpi=720,bbox_inches='tight')
    
    # also plot histograms of the reduced shears at the median-masses
    g1_residual, g2_residual = multipeak_residuals(median_masses,**complete_sample)
    
    # bin in this interval 
    bins = np.linspace(-2,2,101)
    
    fig,ax = pl.subplots()
    ax.set_xlabel(r"$g_{(i)} - g^{(NFW)}_{(i)}$")
    ax.set_ylabel(" Count ")
    
    ax.hist(g1_residual,bins=bins,label=r"$g_1 - g^{(NFW)}_{1}$",histtype='step')
    ax.hist(g2_residual,bins=bins,label=r"$g_2 - g^{(NFW)}_{2}$",histtype='step')
    ax.legend()
    fig.savefig(output_directory + 'shear_residuals.png',dpi=720,bbox_inches='tight')
    
    # get theoretical shears from the residuals
    g1_observed = complete_sample['g1_observed']
    g2_observed = complete_sample['g2_observed']
    
    
    # draw this FOR EACH PEAK, with profiles centered at the location of the peak. From g1/g2 this is actually an easy calculation; this won't work well for multiple large peaks
    for i in range(1): #len(masses[0,:])):
        
        center = complete_sample['cluster_centers'][i]
        delta_x = np.arange(100,3000,10)  # compute the profile out to 3000" ~11.4e3px, out to the edge of the central patches
        delta_y = np.arange(100,3000,10)
        radius_theory = np.sqrt(delta_x**2 + delta_y**2)
        delta_sigma = shear_nfw(median_masses[i]*1e14,delta_x = delta_x*np.pi/( 3600*180 ),delta_y = delta_y*np.pi/( 3600*180 ),cosmo=cosmo,zL = complete_sample['zL'],zS=None)
        delta_sigma_upper = shear_nfw(upper_masses[i]*1e14,delta_x = delta_x*np.pi/( 3600*180 ),delta_y = delta_y*np.pi/( 3600*180 ),cosmo=cosmo,zL = complete_sample['zL'],zS=None)
        delta_sigma_lower = shear_nfw(lower_masses[i]*1e14,delta_x = delta_x*np.pi/( 3600*180 ),delta_y = delta_y*np.pi/( 3600*180 ),cosmo=cosmo,zL = complete_sample['zL'],zS=None)
        
        # compute the E-mode and B-mode shears per-object
        delta_x = (complete_sample['x_pos'] - center[0])*resolution
        delta_y = (complete_sample['y_pos'] - center[1])*resolution
        radius = np.sqrt(delta_x**2 + delta_y**2) 
        theta = np.arctan2(delta_y,delta_x)
        gT = (-g1_observed*np.cos(2*theta) - g2_observed*np.sin(2*theta)) * sigma_crit(zS = complete_sample['zS'],zL=complete_sample['zL'],cosmo=cosmo)
        gX =  (g1_observed*np.sin(2*theta) - g2_observed*np.cos(2*theta)) * sigma_crit(zS = complete_sample['zS'],zL=complete_sample['zL'],cosmo=cosmo)
        
        # bin the observations in radii
        num_bins=18
        stderr = lambda x : np.std(x)/np.sqrt(len(x))
        meanT,bin_edges,binnumber = stats.binned_statistic(radius,gT,bins=num_bins,statistic='median')
        stdT,tmp,binnumber = stats.binned_statistic(radius,gT,bins=bin_edges,statistic=stderr)
        meanX,tmp,binnumber = stats.binned_statistic(radius,gX,bins=bin_edges,statistic='median')
        stdX,tmp,binnumber = stats.binned_statistic(radius,gX,bins=bin_edges,statistic=stderr)
        
        fig,ax = pl.subplots()
        ax2 = ax.twiny()
        #ax.plot(mean_r,mean_gX,'.',label='gX')
        bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in range(len(bin_edges)-1)]
        ax.errorbar(x=bin_centers,y=meanT/1e14,yerr=stdT/1e14,fmt='.',capsize=3,label=r"$\Delta\Sigma_{+}$")
        ax.errorbar(x=bin_centers,y=meanX/1e14,yerr=stdX/1e14,fmt='.',capsize=3,label=r"$\Delta\Sigma_{X}$")
        ax.fill_between(radius_theory,delta_sigma_upper/1e14,delta_sigma_lower/1e14,alpha=0.5,label=r"$\Delta\Sigma_{NFW}$",color='C0')
        ax.hlines(y=0, xmin=0, xmax=4000, color='gray', linestyle='--', alpha=0.5)
        ax.legend()
        
        # set limits, these are hard-coded but generally should work well
        #TODO is there a good data-driven way of assigning these that will look nice?
        ax.set_xlim((200,3050)) # this is generally the edge of our deepCoadd data
        ax.set_ylim(-1.5) # just to cut-off the lower errorbars
        
        # lower-axis in kpc, upper-axis in arcsec
        # arcsec
        tick_locations = (180*3600/np.pi)*np.array([500,1000,1500,2000,2500,3000])/(cosmo.angular_diameter_distance(complete_sample['zL']).value * 1e3)
        # kpc
        tick_labels = np.array([500,1000,1500,2000,2500,3000])
        
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(tick_locations)
        ax2.set_xticklabels(tick_labels)
        ax2.set_xlabel(" r [kpc] ")
        ax.set_xlabel(" r [\"] ")
        ax.set_ylabel(r" $\Delta\Sigma(r)$ $[10^{14} M_\odot Mpc^{-2}]$ ")
        
        fig.savefig(output_directory + 'delta_sigma_profile.png',dpi=720,bbox_inches='tight')

    return

    
if __name__ == '__main__':
    
    # collecting arguments from cln
    if len(sys.argv) == 8:

        table_filename = sys.argv[1]
        Map_table_filename = sys.argv[2]
        cluster_name = sys.argv[3]
        sample_spacing = int(sys.argv[4]) #TODO this and lower-patch can be removed when mass_map outputs either on-sky OR in 00-1111 WCS
        output_directory = sys.argv[5]
        cores = int(sys.argv[6])
        instrument = sys.argv[7]
        lower_patch = (3,3)
        
    elif len(sys.argv) == 9:
    
        table_filename = sys.argv[1]
        Map_table_filename = sys.argv[2]
        cluster_name = sys.argv[3]
        sample_spacing = int(sys.argv[4]) 
        output_directory = sys.argv[5]
        cores = int(sys.argv[6])
        instrument = sys.argv[7]
        lower_patch = sys.argv[8].split(',')
    
    else:
    
        print("python mass_map.py table coadd grid_resolution filter output_directory cores [OPTIONAL: LOWER_PATCH]")
        raise Exception("Improper Usage! Correct usage: python mass_map.py table coadd grid_resolution filter output_directory cores [OPTIONAL: LOWER_PATCH]")
    
    # load and cut the table
    table = ascii.read(table_filename)
    table = load_quality_cuts(table)
    basename = Path(table_filename).stem
    table.write(output_directory + basename + '_mass_fit_cut.csv',format='ascii.csv',overwrite=True)

    peak_table = ascii.read(Map_table_filename)
    
    # by default, run 1000 bootstraps
    realizations = 1000
    
    # load resolution
    resolution_dict = load_resolution()
    resolution = resolution_dict[instrument] # ("/px)
    rad_per_px = (resolution) * np.pi/(180*3600)
    
    #TODO we could automatically select the highest 'unique' N-peaks, something fun for someone to work on; for now I've hard-coded it to 1
    # get px-coordinates of the maximum-SN peak(s)
    num_peaks = 1
    peak_SN = np.max(peak_table['SN_peak'])
    cluster_peak_index = np.argmin(np.abs(peak_SN - peak_table['SN_peak']))
    cluster_centers = [(int(lower_patch[0])*4000 + sample_spacing*(peak_table['x_sn_max'][cluster_peak_index] + 1/2),int(lower_patch[1])*4000 + sample_spacing*(peak_table['y_sn_max'][cluster_peak_index] + 1/2))]

    # we need to know the cluster for the lens-redshift
    ned_result = Ned.query_object(cluster_name)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    zL = ned_result[0]['Redshift']
    
    # preparing everything for the minimization
    x_pos = table['x']
    y_pos = table['y']
    g1_observed = table['g1']
    g2_observed = table['g2']
    zS = table['z_phot'] #z_b in our old catalogs

    #TODO this should probably be moved to the naive shear_calibration
    if 'shape_weight' not in table.colnames:
        weights = np.ones(len(table))/(0.365**2) # shape-weight is an inverse-variance, from previous studies e_rms ~ 0.365 
    else:
        weights = table['shape_weight'].value
    
    complete_sample = {'x_pos':x_pos,'y_pos':y_pos,'g1_observed':g1_observed,'g2_observed':g2_observed,'zS':zS,'zL':zL,'num_peaks':num_peaks,'cluster_centers':cluster_centers,'rad_per_px':rad_per_px}
    
    # each subprocess needs an individual seed, otherwise rng-sequence will repeat in parallel subprocesses
    # defined in main() to use global variables above, returns M200 and the reduced-chi2
    def wrapper(index):
    
        print(index)
        local = np.random.RandomState(index)
        sample = local.randint(low=0,high=len(zS),size=int(len(zS)))
        
        # defining functions inside functions isn't a great thing to-do... but otherwise propagating args through minimize/pool is a pain
        #TODO there is probably a cleaner way of wrapping all this in one method to pass to Pool...
        def chi2(log10_M200):
            g1_res,g2_res = multipeak_residuals([10**power for power in log10_M200],zS[sample],x_pos[sample],y_pos[sample],g1_observed[sample],g2_observed[sample],zL=zL,num_peaks=num_peaks,cluster_centers=cluster_centers,rad_per_px=rad_per_px)
            chi2 = np.sum( ( (g1_res**2) + (g2_res**2) )*weights[sample])
            return chi2
        
        test = minimize(chi2,[14] * num_peaks,bounds=[(12,16)] * num_peaks)
        
        return 10**test.x, test.fun/(len(zS) - num_peaks)

    
    with Pool(cores) as p:
        outputs = p.map(wrapper,range(realizations))

    #TODO from N-peaks, what is the best way to unpack the outputs?
    # right now I'm only unpacking a single-peak, but most of the code is in-place to run statistics on multiple peaks
    mass = []
    chi2 = []
    for i in range(len(outputs)):
        chi2.append(outputs[i][1])
        mass.append(outputs[i][0][0])
    
    draw_plots(np.array(mass),[cluster_name],reduced_chi2=np.array(chi2),complete_sample=complete_sample,output_directory=output_directory,resolution=resolution)


