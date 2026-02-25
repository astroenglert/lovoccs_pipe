
import sys

import numpy as np
import matplotlib.pyplot as pl

import astropy.units as u

from astropy.io import ascii, fits
from astropy.table import Table

from astropy.wcs import WCS
from astropy.visualization import ManualInterval
from astropy.visualization import AsinhStretch

'''

Relatively simple helper script here, it loads the table of mass_peaks and identifies the largest SN-peak (and renders the largest peak as a function of aperture). Then, it takes an input fits file, inverts it, and overlays the contours from mass_map onto it. 

'''

#TODO this function may need to be overhauled as we start searching for peak-SN with trees and things
def find_peak_sn(peak_table,map_output_directory,aperture='schirmer'):
    '''
    
    TODO
    
    '''
    
    peak_index = np.argmax(peak_table['SN_peak'])
    aperture_size = peak_table['aperture_size'][peak_index]
    
    filepath = map_output_directory + 'Map_' + aperture + '_%s.fits'%(int(aperture_size))
    
    return filepath


# our standard asinh stretching
def display(mat, vmin=-0.1, vmax=250):
    interval = ManualInterval(vmin, vmax)
    mat_1 = interval(mat)
    stretch = AsinhStretch(0.0004)
    mat_2 = stretch(mat_1)
    
    return mat_2


if __name__ == '__main__':
    
    if len(sys.argv) == 5:
    
        image_filepath = sys.argv[1]
        peak_table_filepath = sys.argv[2]
        map_output_directory = sys.argv[3]
        output_filename = sys.argv[4]
        xray_filepath = None
        draw_xray = False
    
    elif len(sys.argv) == 6:
        
        image_filepath = sys.argv[1]
        peak_table_filepath = sys.argv[2]
        map_output_directory = sys.argv[3]
        output_filename = sys.argv[4]
        xray_filepath = sys.argv[5]
        draw_xray = True
    
    else:
        raise Exception("Improper usage! Correct usage: python image_filepath peak_table_filepath map_output_directory output_directory [OPTIONAL: xray_filepath]")
    
    # load files from disk
    image_hdul = fits.open(image_filepath)
    peak_table = ascii.read(peak_table_filepath)
    
    # load the peak_table and find the peak aperture
    map_filepath = find_peak_sn(peak_table,map_output_directory)
    contour_hdul = fits.open(map_filepath)
    
    # load data/wcs
    image_data = fits.getdata(image_filepath)
    image_wcs = WCS(image_hdul[0].header)
    
    # load the wcs and compute the SN
    contour_wcs = WCS(contour_hdul[0].header)
    Map_E = contour_hdul[0].data
    Map_V = contour_hdul[2].data
    Map_SN = Map_E/np.sqrt(Map_V)
    
    # now draw the image with the appropriate WCS
    pl.figure(figsize=(8,8))
    ax = pl.subplot(projection=image_wcs,coord_meta={'unit':(u.deg,u.deg)},aspect='equal')
    ax.imshow(1-display(image_data),cmap='gray')
    
    # next overlay the contours
    cs = ax.contour(Map_SN,transform=ax.get_transform(contour_wcs),levels=np.arange(-3,9),alpha=0.8)
    ax.clabel(cs, inline=True, fontsize=11, fmt='%.0f')
    
    # optionally overlay x-ray contours
    if draw_xray:
        
        x_data = fits.getdata(xray_filepath)
        x_wcs = WCS(fits.open(xray_filepath)[0].header)
        
        ax.contour(x_data,transform=ax.get_transform(x_wcs),levels=np.linspace(np.min(x_data),np.max(x_data),10),colors='salmon')
    
    # formatting the axes
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_major_formatter('d.d')
    lat.set_major_formatter('d.d')
    lon.set_axislabel('RA')
    lat.set_axislabel('DEC')
    pl.savefig(output_filename,dpi=720)
    pl.close()
    
    
