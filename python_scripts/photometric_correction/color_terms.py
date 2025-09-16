# color_term correction
# apply color_terms and calculate the bias in residuals, then save that to a file

import sys
import os

from astropy.table import Table

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl

# homebrew modules below
from ..configs.photometric_correction_config import instrument_headers, std_star_headers, colorterm_dictionary, std_star_filename


def get_instrument_headers(catalog_name):
    '''
    Outputs a dictionary of column_names for a given catalog.
    
    Args:
      catalog_name: String; one of decam, des, sm, ps1, sdss from which to load the proper headers for a reference catalog
    
    Returns:
        output_dict: Dictionary; a dictionary which maps ra/dec and magnitudes to the appropriate column of a reference catalog
    
    '''

    if catalog_name in instrument_headers.keys():
        output_dict = instrument_headers[catalog_name]
    else:
        raise Exception(f'Dictionary not found for {catalog_name}')
    
    return output_dict


def get_std_star_sed_mag_headers(catalog_name):
    '''
    Unfortunately, our table of the photometry for reference stars uses a different convention for the column headers, hence creating this function. For now, this is a temporary function which outputs a hand-coded dictionary of column names. eventually we need to replace this with a full solution which loads the appropriate configs from jsons or some other format. 
    
    Args:
      catalog_name: String; one of decam, des, sm, ps1, sdss from which to load the proper headers for a reference catalog
    
    Returns:
        output_dict: Dictionary; a dictionary which maps ra/dec and magnitudes to the appropriate column of a reference catalog
    
    '''

    #output_dict = None
    if catalog_name in std_star_headers.keys():
        output_dict = std_star_headers[catalog_name]
    else:
        raise Exception(f'Dictionary not found for {catalog_name}')
    
    return output_dict

# alright... now we have some data saved to csv's!

# ALGORITHM OUTLINE
# 1. load color_terms and standard stars from disk
# 2. apply color_terms to refcat(to correct for differences in bandpasses)
# 3. compute residuals between color-term corrected catalog and refcat
# 4. assume additive/linear bias and save it to an output file to be applied later

# will depend on refcat instrument!
def load_color_terms(catalog_instr,refcat_instr):
    '''
    Load color_terms for converting between catalog_instr and refcat_instr magnitudes, this will have to be overhauled eventually but for now I'll hard-code a unified dictionary with the locations of the color-terms saved.
    
    Args:
        catalog_instr: string; string specifying the instrument of the catalog
        refcat_instr: string; string specifying the instrument of the refcat
    
    Returns:
        colorterms: Astropy Table; Table with header specifying the band and n-th row specifying n-th order colorterm
    '''
    
    if refcat_instr not in colorterm_dictionary.keys():
        raise Exception(f'No colorterms for {refcat_instr}!')
    else:
        colorterms = Table.read(colorterm_dictionary[refcat_instr])
        
    return colorterms

# will depend on refcat instrument!
def load_std_stars(instr):
    '''
    Load the photometry of standard stars for a given instrument
    
    Args:
        instr: string; string specifying the instrument 
    
    Returns:
        std_star_table: Astropy Table; Table containing the flux of reference stars for the catalog_instr and refcat_instr
    '''
    
    if std_star_filename == None:
        return None
    
    # loading the complete std_star_table and the corresponding headers
    std_star_table_complete = Table.read(std_star_filename)
    std_star_headers_dict = get_std_star_sed_mag_headers(instr)
    
    # collect the bands which need to be written to the table (not None in the dictionary)
    output_me = []
    for band in std_star_headers_dict.keys():
    
        if std_star_headers_dict[band] == None:
            continue
        else:
            output_me.append(band)
    
    # now create a new table which will have updated headers
    std_star_table = Table()
    
    for band in output_me:
        std_star_table[band] = std_star_table_complete[std_star_headers_dict[band]]
    
    return std_star_table

# will (maybe) depend on refcat instrument!
def apply_color_terms(matched_catalog,catalog_instr,refcat_instr,colorterms,catalog_tag='_cat',refcat_tag='_ref',color=['g_psf_mag','i_psf_mag'],single_band=None):
    '''
    Apply color_terms using provided array as coefficients. Applies the correction to all bands where catalog/refcat share a band according to get_instrument_dict()
    
    Args:
        matched_catalog: Astropy Table; Table of matched sources (between a catalog and a refcat)
        catalog_instr: string; string specifying the instrument of the catalog
        refcat_instr: string; string specifying the instrument of the refcat
        colorterms: Astropy Table; Table with header specfiying band and n-th row specifying n-th order colorterm
        color: Array; array specifying the mags to use for color, color[0] - color[1], defaults to g-i
        single_band: string; Optional string for applying color-terms to only one-specific band
        
    Returns:
        residuals: Astropy Table; residuals btwn catalog/refcat
        selected: Numpy Array; array of the sources used for computing bias (aligned with matched_catalog)
    '''
    
    # collecting catalog/refcat headers
    catalog_headers = get_instrument_headers(catalog_instr)
    refcat_headers = get_instrument_headers(refcat_instr)
    
    # select decam bands which are in refcat
    keys = list(refcat_headers.keys())
    mags = []
    for entry in keys:
        # _mag follows keys containing mag, and refcat_headers is None if the band/mag is not in refcat
        if entry[-4:] == '_mag' and refcat_headers[entry] is not None:
            mags.append(entry)
    
    #TODO somehow we need to write this condition in a smarter way
    # right now this is here to prevent a crash when reaching color-selection since legacy lacks i-band
    # for the special-case of des/legacy, we can return zero by default
    # to be fair, we still apply color-cut but to the observed-catalog (since legacy doesn't include i-band as of dr9)
    if (refcat_instr == 'legacy') or (refcat_instr == 'des'):
        
        # for des/legacy, use decam-colors for cuts
        #TODO should we apply the color-cuts here? What are other people's thoughts?
        select = (matched_catalog[catalog_headers[color[0]]+catalog_tag] - matched_catalog[catalog_headers[color[1]]+catalog_tag]) > 0.
        select &= (matched_catalog[catalog_headers[color[0]]+catalog_tag] - matched_catalog[catalog_headers[color[1]]+catalog_tag]) < 0.7
        select_catalog = matched_catalog[select]
        
        # then compute the difference between the two, just w.out color-terms
        residuals = Table()
        for entry in mags:
            band = entry[0]
            residuals[band] = (select_catalog[catalog_headers[entry] + catalog_tag] - select_catalog[refcat_headers[entry] + refcat_tag])
        
        return residuals, select 
    
    # our color_terms are only good for 0 < g-i < 0.7; sufficient to get the linear bias
    select = (matched_catalog[refcat_headers[color[0]]+refcat_tag] - matched_catalog[refcat_headers[color[1]]+refcat_tag]) > 0.
    select &= (matched_catalog[refcat_headers[color[0]]+refcat_tag] - matched_catalog[refcat_headers[color[1]]+refcat_tag]) < 0.7
    select_catalog = matched_catalog[select]
    
    color_zero = select_catalog[refcat_headers[color[0]] + refcat_tag]
    color_one = select_catalog[refcat_headers[color[1]] + refcat_tag]
    
    # if single_band is defined, override mags to contain a single entry
    if single_band is not None:
        mags = [single_band]
    
    # define ct_function and colorterms array, then apply the correction and get the residuals
    residuals = Table()
    for entry in mags:
    
        # band is alawys the first character
        band = entry[0]
        
        # colorterms have headers from band
        coeffs = colorterms[band]
        ct_function=lambda c: coeffs[0] + coeffs[1]*c + coeffs[2]*(c**2) + coeffs[3]*(c**3)
        
        # residuals are (DECam - refcat) - ColorTerm
        residuals[band] = (select_catalog[catalog_headers[entry] + catalog_tag] - select_catalog[refcat_headers[entry] + refcat_tag]) - ct_function(color_zero - color_one)
    
    return residuals, select

def draw_residuals(matched_catalog,catalog_instr,refcat_instr,band,colorterms,catalog_tag='_cat',refcat_tag='_ref',color=['g_psf_mag','i_psf_mag']):
    '''
    Return a figure with the residuals (catalog - refcat) plotted and a curve showing the theoretical color-terms.
    
    Args:
        matched_catalog: Astropy Table; Table of matched sources (between a catalog and a refcat)
        catalog_instr: string; string specifying the instrument of the catalog
        refcat_instr: string; string specifying the instrument of the refcat
        band: string; string specifying the band to plot/correct
        colorterms: Astropy Table; Table with header specfiying band and n-th row specifying n-th order colorterm
        color: Array; array specifying the mags to use for color, color[0] - color[1], defaults to g-i
        
    Returns:
        ax: Matplotlib Axes; A pair of axes containing a plot of the residuals v. color
    
    '''
    
    fig,ax = pl.subplots()
    
    # collecting catalog/refcat headers
    catalog_headers = get_instrument_headers(catalog_instr)
    refcat_headers = get_instrument_headers(refcat_instr)
    filter_name = band[0]
    ct0_filter_name = color[0][0]
    ct1_filter_name = color[1][0]
    
    # first load in all data points
    ct_band0_refcat = matched_catalog[refcat_headers[color[0]] + refcat_tag]
    ct_band1_refcat = matched_catalog[refcat_headers[color[1]] + refcat_tag]
    band_refcat = matched_catalog[refcat_headers[band] + refcat_tag]
    band_catalog = matched_catalog[catalog_headers[band] + catalog_tag]
    
    # filter-out NaN/inf
    select_plot = np.isfinite(ct_band0_refcat)
    select_plot &= np.isfinite(ct_band1_refcat)
    select_plot &= np.isfinite(band_refcat)
    select_plot &= np.isfinite(band_catalog)
    
    # draw the observations
    im = ax.scatter((ct_band0_refcat - ct_band1_refcat)[select_plot],
                    (band_catalog - band_refcat)[select_plot],
                    c=band_refcat[select_plot],
                    s=6.,
                    marker='o',
                    cmap='viridis',
                    alpha=0.8,
                    label='observation',
                    )
    
    # next load in the standard stars
    refcat_std_stars = load_std_stars(refcat_instr)
    catalog_std_stars = load_std_stars(catalog_instr)
    
    # check that the reference stars exist, if-so draw them as well
    if (refcat_std_stars is not None) and (catalog_std_stars is not None):
        refcat_std_ct0 = refcat_std_stars[color[0]]
        refcat_std_ct1 = refcat_std_stars[color[1]]
        refcat_std_band = refcat_std_stars[band]
        catalog_std_band = catalog_std_stars[band]
    
        # draw the standard stars
        ax.scatter((refcat_std_ct0 - refcat_std_ct1),
                   (catalog_std_band - refcat_std_band),
                   c='k',
                   alpha=0.5,
                   marker='^',
                   label='model',
                   )
    
    # draw the theoretical color-terms    
    coeffs = colorterms[filter_name]
    ct_function=lambda c: coeffs[0] + coeffs[1]*c + coeffs[2]*(c**2) + coeffs[3]*(c**3)
    theory_x = np.linspace(0,0.7,100)
    theory_y = ct_function(theory_x)
    
    ax.plot(theory_x,theory_y, 'k--')
    
    fig.colorbar(im)
    ax.set_xlim([-0,0.7])
    ax.set_ylim([-0.2,0.2])
    ax.set_ylabel(f'{filter_name}_{catalog_instr} - {filter_name}_{refcat_instr}')
    ax.set_xlabel(f'{ct0_filter_name}_{refcat_instr} - {ct1_filter_name}_{refcat_instr}')
    ax.legend('lower right')
    
    # draw color_term residuals
    return fig,ax

if __name__ == '__main__':
    
    if len(sys.argv)==5:
    
        # collecting arguments from cln
        matched_catalog_filename = sys.argv[1]
        catalog_instr = sys.argv[2] # usually decam
        refcat_instr = sys.argv[3] # ps1, sm, sdss, des
        residual_filename = sys.argv[4]
        single_band = None # if this isn't specified, run color_terms on all bands shared with the refcat
        
    elif len(sys.argv) == 6:
    
        # collecting arguments from cln
        matched_catalog_filename = sys.argv[1]
        catalog_instr = sys.argv[2] # usually decam
        refcat_instr = sys.argv[3] # ps1, sm, sdss, des
        residual_filename = sys.argv[4]
        single_band = sys.argv[5] # optional argument for running on a single-band (usually u-band)
        
    else:
        print("python color_terms.py matched_catalog catalog_instr refcat_instr residual_filename [OPTIONAL: band]")
        raise Exception("Improper Usage! Correct usage: python color_terms.py matched_catalog catalog_instr refcat_instr residual_filename [OPTIONAL: band]")
    
    # if residual_filename exists, load it, otherwise create a new table
    if os.path.exists(residual_filename):
        residual_table = Table.read(residual_filename,format='ascii.csv')
    else:
        residual_table = Table()
    
    # collect dictionaries for headers
    refcat_headers = get_instrument_headers(refcat_instr)
    catalog_headers = get_instrument_headers(catalog_instr)
    
    # load a matched_catalog and the appropriate colorterms
    matched_catalog = Table.read(matched_catalog_filename,format='ascii.csv')
    colorterms = load_color_terms(catalog_instr,refcat_instr)
    
    # apply the color-terms and collect the residuals for all bands btwn catalog/refcat
    residuals,selected = apply_color_terms(matched_catalog,catalog_instr,refcat_instr,colorterms,single_band=single_band)
        
    filters = residuals.colnames
    
    for i in range(len(filters)):
    
        band = filters[i] + '_psf_mag' # assume that we're using psf_mag for correction by default
        
        # draw and save the residual figure
        fig,ax = draw_residuals(matched_catalog,catalog_instr,refcat_instr,band,colorterms)
        fig.savefig("%s_%s_%s.png"%(os.path.splitext(residual_filename)[0], band,refcat_instr),bbox_inches='tight')
        
        # now let's collect the bias using the residuals, apply the magnitude-cut here
        refcat_tag = '_ref'
        refcat_selected_mag = matched_catalog[refcat_headers[band] + refcat_tag][selected]
        
        if filters[i] == 'u':
            magnitude_cut = refcat_selected_mag > 14
        else:
            magnitude_cut = refcat_selected_mag > 16
        total_selected = np.sum(selected)
        print(f'Using {total_selected} to compute a linear-bias in the residuals')
        # now collect the final residuals
        final_residuals = residuals[filters[i]][magnitude_cut]
        
        # save a histogram of the residuals
        fig,ax = pl.subplots()
        n,b,_ = ax.hist(final_residuals, bins=np.linspace(-0.15, 0.15, 31)) 
        mid = (b[1:]+b[:-1])/2.
        print("peak: ", mid[np.argmax(n)])

        fig.savefig("%s_%s_bias_hist.png"%(os.path.splitext(residual_filename)[0], filters[i]),bbox_inches='tight' )

        # update the residuals table
        median_bias = np.nanmedian(final_residuals)
        std_bias = np.nanstd(final_residuals)/np.sqrt(np.sum(np.isfinite(final_residuals)))
        
        print(f'Median Bias: {median_bias:.3f}, Std Bias: {std_bias:.3f}')
        residual_table[filters[i]] = [median_bias,std_bias]
    
    residual_table.write(residual_filename, format="ascii.csv", overwrite=True)
    
    
    
    
    
    
        


