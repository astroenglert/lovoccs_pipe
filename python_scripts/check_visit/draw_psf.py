import os
import time
import math
import sys
from multiprocessing import Pool
import matplotlib.pyplot as pl

# disable implicit threading so tasks can run in parallel without fighting for threads
os.environ['MKL_NUM_THREADS'] = '1'

import numpy as np
import astropy as ap
import pandas as pd

from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord

from lsst.daf.butler import Butler

# helper function to create colorbars without adjusting aspect-ratio
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

def draw_psf(visit,band='r',fwhm_cut=3.8,ellip_cut=0.1,res=120):
    '''
    Draw a model of the PSF across a visit
    
    Args:
      visit: int; visit number to draw
      band: string; band for this visit
      fwhm_cut: float; maximum FWHM to pass QC's
      ellip_cut: float; maximum ellip to pass QC's
      res: float; resolutiono f the final figure
    
    Returns:
      None
    
    '''

    visit_summary_filepath = "check_visit/summary_tables/{band}_visit_{exp}_summary.csv".format(band=band,exp=visit)
    
    visit_wcs_filepath = "check_visit/summary_tables/{band}_visit_{exp}_wcs.csv".format(band=band,exp=visit)
    
    # loading the output from check_visit
    visit_summary = pd.read_csv(visit_summary_filepath,index_col=0)
    
    # loading wcs-coords
    visit_wcs = pd.read_csv(visit_wcs_filepath,index_col=0)
    
    # removing the bad-detectors, at least 2,31,61 (need to ask about 62)
    visit_bad_detectors = (visit_summary['detector'] == 2) | (visit_summary['detector'] == 31) |  (visit_summary['detector'] == 61)
    visit_summary = visit_summary[np.invert(visit_bad_detectors)]
    
    wcs_bad_detectors = (visit_wcs.index == 2) | (visit_wcs.index == 31) |  (visit_wcs.index == 61)
    visit_wcs = visit_wcs[np.invert(wcs_bad_detectors)]
    ra_dec_pairs = visit_wcs[['ra','dec']].values
    
    # compute per-detector statistics
    detectors = np.unique(visit_summary['detector']).astype(int)
    fwhm_median = []
    ellip_median = []
    
    for det in detectors:
        
        # per-detector median fwhm
        fwhm_median.append( np.median(visit_summary['fwhm'][visit_summary['detector'] == det]) )

        # per-detector median ellip
        e1 = visit_summary['e1'][visit_summary['detector']==det]
        e2 = visit_summary['e2'][visit_summary['detector']==det]
        ellip_median.append( np.median( np.sqrt( e1**2 + e2**2 ) ) )

    # now loading everything else I need from the table
    ra = visit_summary['ra']
    dec = visit_summary['dec']
    fwhm = visit_summary['fwhm']
    e1 = visit_summary['e1']
    e2 = visit_summary['e2']
    ellip = np.sqrt(e1**2 + e2**2)
    u = visit_summary['u']
    v = visit_summary['v']
    
    # first fwhm
    fig,ax = pl.subplots()
    im = ax.scatter(ra,dec,marker='.',c=fwhm,cmap='jet',alpha=0.4)
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("DEC (deg)")
    colorbar(im)
    ax.set_title("Visit {exp} Moffat FWHM".format(exp=visit))
    
    # add per-detector statistics
    for i in range(len(detectors)):
        
        if fwhm_median[i] > fwhm_cut:
            color="Red"
        else:
            color="Black"
            
        ax.text(ra_dec_pairs[i][0],ra_dec_pairs[i][1],'CCD%02d\n%.2f'%(int(detectors[i]), fwhm_median[i]),fontsize='x-small',horizontalalignment='center',verticalalignment='center',color=color)
    
    pl.savefig("check_visit/plots/{band}_visit_{exp}_fwhm.png".format(band=band,exp=visit),bbox_inches='tight',dpi=res)
    pl.close()
    
    # second ellip
    fig,ax = pl.subplots()
    im = ax.scatter(ra,dec,marker='.',c=ellip,cmap='jet',alpha=0.4)
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("DEC (deg)")
    colorbar(im)
    ax.set_title("Visit {exp} Ellipticity".format(exp=visit))
    
    # add per-detector statistics
    for i in range(len(detectors)):
        
        if ellip_median[i] > ellip_cut:
            color="Red"
        else:
            color="Black"
        
        ax.text(ra_dec_pairs[i][0],ra_dec_pairs[i][1],"CCD%02d\n%.2f"%(int(detectors[i]), ellip_median[i]),fontsize='x-small',horizontalalignment='center',verticalalignment='center',color=color)
    
    pl.savefig("check_visit/plots/{band}_visit_{exp}_ellip.png".format(band=band,exp=visit),bbox_inches='tight',dpi=res)
    pl.close()
    
    # third stokes parameters
    fig,ax = pl.subplots()
    im = ax.quiver(ra,dec,u,v,ellip,cmap='jet',alpha=0.4)
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("DEC (deg)")
    colorbar(im)
    ax.set_title("Visit {exp} Ellipticity".format(exp=visit))
    
    # add per-detector statistics
    for i in range(len(detectors)):
        
        # color the text red if this detector is going to fail the corresponding cut
        if ellip_median[i] > ellip_cut:
            color="Red"
        else:
            color="Black"
        
        ax.text(ra_dec_pairs[i][0],ra_dec_pairs[i][1],"CCD%02d\n%.2f"%(int(detectors[i]), ellip_median[i]),fontsize='x-small',horizontalalignment='center',verticalalignment='center',color=color)
    
    pl.savefig("check_visit/plots/{band}_visit_{exp}_ellip_quiver.png".format(band=band,exp=visit),bbox_inches='tight',dpi=res)
    pl.close()

if __name__ == '__main__':
    
    # checking for correct usage
    if len(sys.argv)!=5:
        print('Usage: python this_script.py band cores fwhm_cut ellip_cut')
        sys.exit(1)
    
    band = sys.argv[1]
    cores = int(sys.argv[2])
    fwhm_cut = float(sys.argv[3])
    ellip_cut = float(sys.argv[4])
    
    # create a butler to read/write data with
    butler = Butler('repo/repo',writeable=False)
    
    # load all the visits for a given band (re-using code from check_visit)
    exposure_generator = butler.registry.queryDimensionRecords('exposure',datasets='calexp',collections='DECam/processing/calexp_{band}'.format(band=band))
    exposure_ids = []
    
    for exposure in exposure_generator:
        exposure_ids.append(exposure.id)
    
    exposure_ids = np.unique(np.array(exposure_ids))
     
    # split into workers and start plotting visits (again, re-using code from check_visit)
    num_processes = cores
    groups = math.ceil(len(exposure_ids)/num_processes)
    for i in range(groups):
        
        print("Now processing group {current} of {total}".format(current=i+1,total=groups))
        if i == groups-1:
            process_visits = exposure_ids[i*num_processes:]
        else:
            process_visits = exposure_ids[i*num_processes:(i+1)*num_processes]
        
        # collect the correct arguments
        arguments = []
        for visit in process_visits:
            arguments.append((visit,band,fwhm_cut,ellip_cut))
            
        # create a pool for parallelization
        with Pool(len(process_visits)) as pool:
            output = pool.starmap(draw_psf, arguments)

