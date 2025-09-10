# helper function to create either a fits or png for the given patch range
import os
import time
import math
import glob
import sys

import numpy as np

import astropy as ap
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.visualization import ManualInterval
from astropy.visualization import AsinhStretch
from astropy.visualization import LogStretch
from astropy.visualization import make_lupton_rgb
from astropy.table import Table

import pandas as pd
import matplotlib.pyplot as pl

from lsst.daf.butler import Butler

# by default we clip and stretch the data before displaying it
def display(mat, vmin, vmax, tag):
    interval = ManualInterval(vmin, vmax)
    mat_1 = interval(mat)
    stretch = AsinhStretch(0.0004)
    mat_2 = stretch(mat_1)

    return mat_2

# Let's also use two stretches, one for color-correct (lupton) and another for not color-correct (but nice looking)
def draw_rgb(R,G,B,out,tag='deepCoadd',lupton=False):
    '''
    Helper function to draw an RGB image of the data
    
    Args:
      R/G/B: string; string specifying the band to load into each channel
      patches: string; the range of patches (e.g. '33-88')
      tag: string; tag to append to the end of the output file
      lupton: bool; enable Lupton coloring?
    
    Returns
      None
    
    '''
    
    print("Loading FITS image...")
    if os.path.exists(R) and os.path.exists(G) and os.path.exists(B):
        print("Files all exist!")
    else:
        print("Some files could be missing...\nExiting...")
        return False
    data_1 = (pyfits.getdata(R, 0))
    data_2 = (pyfits.getdata(G, 0))
    data_3 = (pyfits.getdata(B, 0))
    
    # for pretty-pictures, -0.1 and 250 seem to work best
    # but for science -0.06 and 190 best showcase our issues with bckg-subtr and colors
    
    # if Lupton coloring is enabled, we use a different transformation
    print("Scaling image...")
    
    if lupton:
        data_all = make_lupton_rgb(data_1,data_2,data_3,minimum=-0.1,Q=8,stretch=1)
        print("Writing lupton_%s..."%out)
        pl.imsave("lupton_%s"%out, 
            np.flipud(data_all),
           )

    else:
        data_1_scaled = display(data_1, -0.1, 250, 'fig1')*255
        data_2_scaled = display(data_2, -0.1, 250, 'fig2')*255
        data_3_scaled = display(data_3, -0.1, 250, 'fig3')*255
        
        data_1_scaled = np.flipud(data_1_scaled)
        data_2_scaled = np.flipud(data_2_scaled)
        data_3_scaled = np.flipud(data_3_scaled)
        data_all = np.dstack((
                                data_1_scaled, 
                                data_2_scaled, 
                                data_3_scaled,
                            )).astype(np.uint8)
        print("Writing %s..."%out)
        pl.imsave(out, 
            data_all,
           )
    return True


if __name__=='__main__':
    
    # one use case, provide filepaths to the three fits you want to combine and a binary flag for Lupton (optional, 1 enables 0 disables)  
    if not (len(sys.argv)==6 or len(sys.argv)==5):  
        print("Improper usage, try either: ")
        print("python draw_patches_rgb.py R G B filename OR python draw_patches_rgb.py R G B filename Lupton")
        sys.exit(1)
    
    R_path = sys.argv[1]
    G_path = sys.argv[2]
    B_path = sys.argv[3]
    filename = sys.argv[4]
    
    if len(sys.argv)==6:
        Lupton = bool(int(sys.argv[5]))
    else:
        Lupton = False
    
    draw_rgb(R_path,G_path,B_path,filename,lupton=Lupton)
    
