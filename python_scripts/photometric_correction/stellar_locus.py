# stellar-locus correction
# plot of u-g v. g-r should be linear, use this to further improve the zero-point for u-band

import sys
import os

import numpy as np

import matplotlib.pyplot as pl

from astropy.table import Table
from scipy.optimize import curve_fit

#homebrew modules here, these will likely change in the future
from .color_terms import get_instrument_headers, get_std_star_sed_mag_headers, load_std_stars

def draw_color_plots(star_catalog,instr,x_color,y_color,ax,kwargs,select=None):
    '''
    Draw a plot of x_color v y_color
    
    Args:
        star_catalog: Astropy Table; Table containing only stars
        instr: String; String specifying the instrument for the star_catalog
        x_color: array; An array specifying the color to be plotted on x, e.g. ['g_psf_mag','r_psf_mag']
        y_color: array; An array specifying the bands to be plotted on u, e.g. ['u_psf_mag','g_psf_mag']
        ax: Matplotlib Axes; Axes to draw the plot on
        kwargs: dictionary; Dictionary containing plot parameters
        select: array; An array of booleans specifying which stars to plot
    
    Returns:
        fig: Matplotlib Figure; Figure containing the plot
        ax: Matplotlib Axes; Axes for the plot
    '''
    
    # loading filter headers
    header_dict = get_instrument_headers(instr)
    
    # if select is None, use all stars in the catalog
    if select is None:
        select = np.ones(len(star_catalog)).astype(bool)
    
    # selecting stars
    select_catalog = star_catalog[select]
    
    # loading magnitudes and computing colors
    x_color_0 = select_catalog[header_dict[x_color[0]]]
    x_color_1 = select_catalog[header_dict[x_color[1]]]
    y_color_0 = select_catalog[header_dict[y_color[0]]]
    y_color_1 = select_catalog[header_dict[y_color[1]]]
    
    x_axis = x_color_0 - x_color_1
    y_axis = y_color_0 - y_color_1
    
    im = ax.scatter(x_axis,y_axis,**kwargs)
    
    return ax,im
    
if __name__ == '__main__':
    
    # check cln arguments, load the catalog, apply cuts as necessary, plot and fit
    if len(sys.argv)!=4:
        print("python this.py star_catalog instr residual_filename")
        raise Exception("Improper usage! Correct usage is python this.py star_catalog catalog_instr residual_filename")
    
    # collecting arguments from cln
    star_catalog_filename = sys.argv[1]
    instr = sys.argv[2] # usually decam, in principle with an updated residual_table other catalogs are also supported
    residual_filename = sys.argv[3]
    
    # collect dictionaries for headers
    headers = get_instrument_headers(instr)
    coltag = '_psf_mag'
    
    # load a matched_catalog and the appropriate colorterms
    star_catalog = Table.read(star_catalog_filename,format='ascii.csv')
    
    # running some preliminary cuts
    star_catalog = star_catalog[star_catalog['extendedness'] < 0.5] # just in-case a non-star-catalog is passed
    star_catalog = star_catalog[star_catalog[headers['g' + coltag]] < 18] # only use sufficiently bright stars
    
    # load the table containing the zp-correction
    residual_table = Table.read(residual_filename,format='ascii.csv')
    file_start = os.path.splitext(residual_filename)[0] # start of output filenames from this script
    
    # load the std-stars
    std_star_table = load_std_stars(instr)
    
    # apply the zp-correction to _psf_mag if the residual exists
    for band in residual_table.colnames:
        # skip any nan/inf/None
        if ~np.isfinite(residual_table[band][0]):
            continue
        
        star_catalog[band + coltag] = star_catalog[band + coltag] - residual_table[band][0]
    
    # if u-band is all nan/inf/None, exit the script since there is nothing to correct!
    if np.sum(np.isfinite(star_catalog['u' + coltag])) == 0:
        
        print("WARNING: u-band is all nan/inf/None, assigning NaN for zp-correction in u and exiting!")
        
        residual_table['u'] = [np.nan, np.nan]
        residual_table.write(residual_filename,format='ascii.csv',overwrite=True)
        sys.exit()
    
    # now draw the plots we need
    #TODO not great coding done here, we can probably wrap each block here into single function-calls
    
    #FIRST: g-r v r-i
    cbar_mag = 'g_psf_mag'
    color_x_0 = 'g_psf_mag'
    color_x_1 = 'r_psf_mag'
    color_y_0 = 'r_psf_mag'
    color_y_1 = 'i_psf_mag'
    band_x0 = color_x_0[0]
    band_x1 = color_x_1[0]
    band_y0 = color_y_0[0]
    band_y1 = color_y_1[0]
    
    fig,ax = pl.subplots()
    
    fmt = {'s':6,'marker':'o','alpha':0.8,'cmap':'viridis','c':star_catalog[headers[cbar_mag]]}
    ax,im = draw_color_plots(star_catalog,instr,[color_x_0,color_x_1],[color_y_0,color_y_1],ax,fmt)
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel(f'{cbar_mag}')
    ax.set_xlabel(f'{band_x0} - {band_x1}')
    ax.set_ylabel(f'{band_y0} - {band_y1}')
    
    # I'll manually draw-on the reference stars since it's just a single line of code
    ax.scatter(std_star_table[color_x_0] - std_star_table[color_x_1],
               std_star_table[color_y_0] - std_star_table[color_y_1],
               alpha=0.5,
               marker='^',
               c='k',
               label='model',
              )
    
    ax.set_xlim((-0.5,0.8))
    ax.set_ylim((-0.3,0.3))
    ax.legend(['Observed Stars','Model Stars'])
    
    fig.savefig(f'{file_start}_{band_x0}-{band_x1}_{band_y0}-{band_y1}',bbox_inches='tight')

    #SECOND: i-z v r-i
    cbar_mag = 'g_psf_mag'
    color_x_0 = 'i_psf_mag'
    color_x_1 = 'z_psf_mag'
    color_y_0 = 'r_psf_mag'
    color_y_1 = 'i_psf_mag'
    band_x0 = color_x_0[0]
    band_x1 = color_x_1[0]
    band_y0 = color_y_0[0]
    band_y1 = color_y_1[0]
    
    fig,ax = pl.subplots()
    
    fmt = {'s':6,'marker':'o','alpha':0.8,'cmap':'viridis','c':star_catalog[headers[cbar_mag]]}
    ax,im = draw_color_plots(star_catalog,instr,[color_x_0,color_x_1],[color_y_0,color_y_1],ax,fmt)
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel(f'{cbar_mag}')
    ax.set_xlabel(f'{band_x0} - {band_x1}')
    ax.set_ylabel(f'{band_y0} - {band_y1}')
    
    # I'll manually draw-on the reference stars since it's just a single line of code
    ax.scatter(std_star_table[color_x_0] - std_star_table[color_x_1],
               std_star_table[color_y_0] - std_star_table[color_y_1],
               alpha=0.5,
               marker='^',
               c='k',
               label='model',
              )
    
    ax.set_xlim((-0.3,0.3))
    ax.set_ylim((-0.2,0.3))
    ax.legend(['Observed Stars','Model Stars'])
    
    fig.savefig(f'{file_start}_{band_x0}-{band_x1}_{band_y0}-{band_y1}',bbox_inches='tight')
    
    #THIRD: u-band correction and g-r v. u-g
    
    # loading colnames and such
    cbar_mag = 'g_psf_mag'
    color_x_0 = 'g_psf_mag'
    color_x_1 = 'r_psf_mag'
    color_y_0 = 'u_psf_mag'
    color_y_1 = 'g_psf_mag'
    band_x0 = color_x_0[0]
    band_x1 = color_x_1[0]
    band_y0 = color_y_0[0]
    band_y1 = color_y_1[0]
    
    # we fit to model-stars in a small region of the color-space
    model_select = (std_star_table[color_y_0] - std_star_table[color_y_1]) > 0.81
    xdata = (std_star_table[color_x_0] - std_star_table[color_x_1])[model_select]
    ydata = (std_star_table[color_y_0] - std_star_table[color_y_1])[model_select]
    func = lambda x, a, b: x*a + b
    popt, pcov = curve_fit(func, xdata, ydata)
    print("Best-fit model parameters:")
    print(popt)
    
    # now draw plots of the fit
    fig,ax = pl.subplots()
    
    # first draw the curve-fit
    theory_x = np.linspace(np.min(xdata),np.max(xdata))
    ax.plot(theory_x,func(theory_x, *popt),'k--',label=['Linear Fit'])
    
    # now the model/observed data-points
    select_plot = np.isfinite( star_catalog[headers[color_x_0]] )
    select_plot &= np.isfinite( star_catalog[headers[color_x_1]] )
    select_plot &= np.isfinite( star_catalog[headers[color_y_0]] )
    select_plot &= np.isfinite( star_catalog[headers[color_y_1]] )
    
    fmt = {'label':['observed'],'s':6,'marker':'o','alpha':0.8,'cmap':'viridis','c':star_catalog[headers[cbar_mag]][select_plot]}
    ax,im = draw_color_plots(star_catalog,instr,[color_x_0,color_x_1],[color_y_0,color_y_1],ax,fmt,select=select_plot)
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel(f'{cbar_mag}')
    ax.set_xlabel(f'{band_x0} - {band_x1}')
    ax.set_ylabel(f'{band_y0} - {band_y1}')
    
    # I'll manually draw-on the reference stars since it's just a single line of code
    ax.scatter(std_star_table[color_x_0] - std_star_table[color_x_1],
               std_star_table[color_y_0] - std_star_table[color_y_1],
               alpha=0.5,
               marker='^',
               c='k',
              )
    
    ax.set_xlim((-0,0.8))
    ax.set_ylim((0.3,2))
    ax.legend(['Linear Fit','Observed Stars','Model Stars'])
    
    filename_start = os.path.splitext(residual_filename)[0]
    fig.savefig(f'{file_start}_{band_x0}-{band_x1}_{band_y0}-{band_y1}',bbox_inches='tight')

    #LASTLY, for computing the linear bias we select from a narrow-region of the observed color-space
    select_bias = ( star_catalog[headers[color_y_0]] - star_catalog[headers[color_y_1]] ) > 0.81
    select_bias &= ( star_catalog[headers[color_y_0]] - star_catalog[headers[color_y_1]] ) < 1.6

    obs_xdata = ( star_catalog[headers[color_x_0]] - star_catalog[headers[color_x_1]] )[select_bias]
    obs_ydata = ( star_catalog[headers[color_y_0]] - star_catalog[headers[color_y_1]] )[select_bias]
    bias = obs_ydata - func(obs_xdata, *popt)
    
    # update the residuals table
    median_bias = np.nanmedian(bias)
    std_bias = np.nanstd(bias)/np.sqrt(np.sum(np.isfinite(bias)))
        
    print(f'Median Bias: {median_bias:.3f}, Std Bias: {std_bias:.3f}')
    residual_table = Table()
    residual_table['u'] = [median_bias,std_bias]
    
    residual_table.write(file_start + '_stellar_locus.csv', format="ascii.csv", overwrite=True)



