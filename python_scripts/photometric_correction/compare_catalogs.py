# comparison plot
# draws plots and computes statistics comparing a matched-catalog

import sys
import os

from astropy.table import Table, join

import numpy as np

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as pl

from scipy.stats import norm

# homebrew modules here
from .color_terms import get_instrument_headers, apply_color_terms, load_color_terms

# for this, I want two use-cases:
# 1. Comparing a single band in a matched catalog
# 2. Comparing all shared bands between the matched catalogs
# draw a set of subplots, one displaying residuals across the field and another the histogram of those

#NOTE: the old version of this script, compare_mag_v3b, centered the residuals at zero by subtracting a median, this doesn't change the stdv but I think it's misleading (since you're no longer rendering the residuals in that case!) so I've removed it from this script (this is line 183 of compare_mag_v3b.py)
if __name__ == '__main__':
    
    if len(sys.argv)==4:
    
        # collecting arguments from cln
        matched_catalog_filename = sys.argv[1]
        catalog_instr = sys.argv[2] # usually decam
        refcat_instr = sys.argv[3] # ps1, sm, sdss, des
        single_band = None # if this isn't specified, run on all bands shared with the refcat
        
    elif len(sys.argv)==5:
    
        # collecting arguments from cln
        matched_catalog_filename = sys.argv[1]
        catalog_instr = sys.argv[2] # usually decam
        refcat_instr = sys.argv[3] # ps1, sm, sdss, des
        single_band = sys.argv[4] # optional argument for running on a single-band (usually u-band), e.g. 'u_psf_mag'
        
    else:
        print("python compare_catalogs.py matched_catalog catalog_instr refcat_instr [OPTIONAL: band]")
        raise Exception("Improper Usage! Correct usage: python compare_catalogs.py matched_catalog catalog_instr refcat_instr [OPTIONAL: band]")
    
    # collect instrument headers
    refcat_headers = get_instrument_headers(refcat_instr)
    catalog_headers = get_instrument_headers(catalog_instr)
    
    # load the matched catalog
    matched_catalog = Table.read(matched_catalog_filename,format='ascii.csv')
    
    # load the color-terms
    colorterms = load_color_terms(catalog_instr,refcat_instr)
    
    # apply_color_terms will automatically pick-out all matched bands or just run on a single-band (and includes the color-cut by default)
    residuals,selected = apply_color_terms(matched_catalog,catalog_instr,refcat_instr,colorterms,single_band=single_band)
    
    # ranges for magnitudes
    mag_ranges = {'g':[15,21],'r':[15,21],'i':[15,21],'z':[15,21],'u':[13,21],'Y':[13,20]}
    
    # Generate dict to write differences to
    band_diff_dict = {}

    # now iterate through the headers which exist in residuals and save the results!
    for band in residuals.colnames:
        
        magmin = mag_ranges[band][0]
        magmax = mag_ranges[band][1]
        
        diffs = residuals[band]
        select_catalog = matched_catalog[selected]
        
        # header tags
        magnitude_type = '_psf_mag'
        ref_tag = '_ref'
        cat_tag = '_cat'
        
        # magnitude cuts (based on the refcat)
        final_cut = (select_catalog[refcat_headers[band + magnitude_type] + ref_tag] < magmax)
        final_cut &= (select_catalog[refcat_headers[band + magnitude_type] + ref_tag] > magmin)
        
        # and one last cut to make sure everything is finite
        final_cut &= np.isfinite(diffs)
        
        # run a check for a zero-array and print a warning if-so
        if np.sum(final_cut) == 0 or ~np.isfinite(np.sum(final_cut)):
            print(f'Catalog for {band}-band is empty after cuts, so no comparisons will be plotted!')
            break
        
        # now collect coordinates/differences with the final cut
        ra_cut = select_catalog[catalog_headers['ra_name'] + cat_tag][final_cut]
        dec_cut = select_catalog[catalog_headers['dec_name'] + cat_tag][final_cut]
        diff_cut = diffs[final_cut]
        id_cut = select_catalog['ID' + cat_tag][final_cut]

        # add table to the dictionary
        band_table = Table({'ID': id_cut, 'ra': ra_cut, 'dec': dec_cut, f'{band}_diff': diff_cut})
        band_diff_dict[band] = band_table
        
        # now draw the pretty-pictures :D
        # I could wrap this in a function... but since this is the only big-step in the script I'll cheat for now
        # all of this is basically copied from Shenming's code
        fig = pl.figure(figsize=(12, 6))
        gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1.5, 1])
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        axs = [ax0, ax1]
        
        # rendering residuals across ra/dec
        fmt = {'c':diff_cut,'s':4,'marker':'o','cmap':'Spectral','vmin':-0.05,'vmax':0.05,'alpha':0.5,'edgecolors':None}
        im = axs[0].scatter(ra_cut,dec_cut,**fmt)
        pl.colorbar(im,ax=axs[0])
        
        # creating the histogram
        bins = np.linspace(-0.25, 0.25, 50*2+1)
        mid_points = (bins[1:] + bins[:-1])/2.
        n, b, p = axs[1].hist(diff_cut, bins=bins, histtype="step", log=True, label=band)
        fig.suptitle("dm - ext: %d stars %.1f<%s<%.1f std %.4f"%(len(diff_cut), magmin, band, magmax, np.std(diff_cut) ) )
        print(f'max, min of diff for {band}: ', np.max(diff_cut), np.min(diff_cut))
        
        #axs[1].plot(mid_points, norm.pdf(mid_points, mean_new, std_new)*(np.sqrt(2.*np.pi)*std_new)*amp_new)
        axs[1].set_ylim([8e-1, 1e4])
        axs[1].set_xlabel("magnitude difference (%s - %s)"%(catalog_instr, refcat_instr) )
        axs[1].set_ylabel("count")
        
        pl.savefig("%s_%s_%s_cmp.png"%(os.path.splitext(matched_catalog_filename)[0], band, refcat_instr) )

    # combine individual band tables into one combined table
    bands_list = list(band_diff_dict.keys())
    combined_table = band_diff_dict[bands_list[0]]

    for band in bands_list[1:]:
        combined_table = join(combined_table, band_diff_dict[band], keys=['ID','ra','dec'], join_type='outer')

    combined_table.write(f'photometric_correction_output/{refcat_instr}_mag_diffs.csv', format="ascii.csv", overwrite=True)
        
        
    
