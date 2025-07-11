import os
import time
import math
import sys
import glob
from multiprocessing import Pool

# disable implicit threading so tasks can run in parallel fighting for threads
os.environ['MKL_NUM_THREADS'] = '1'

import numpy as np
import astropy as ap
import pandas as pd

from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.visualization import ZScaleInterval
from astropy import wcs
from astropy.coordinates import SkyCoord

from lsst.daf.butler import Butler

# create a summary of a given visit
def summary_from_stars(visit,detectors,psf_flag='calib_psf_used',band='r'):
    s = time.time()
    print("Now summarizing exposure number {visit}".format(visit=visit))
    # columns for final dataframe returned at the end
    # not SUPER-optimized or anything... but we don't know the no of detections beforehand
    visit_ra = []
    visit_dec = []
    visit_e1 = []
    visit_e2 = []
    visit_u = []
    visit_v=[]
    visit_fwhm = []
    visit_det = []
    psf_candidates_index = []

    ra_dec_pairs = []
    
    # for det in detectors:
    for det in detectors:
    
        # first get the hdul
        hdul = fits.open(butler.getURI('calexp',dataId={'instrument' : 'DECam', 'band':band, 'detector':det, 'visit':visit},collections='DECam/processing/calexp_{band}'.format(band=band)).ospath)
            
        # fits images are a little easier to work with than the internal lsst.afw.image objects
        calexp = hdul[1].data
   
        # also loading wcs/detector
        w = wcs.WCS(hdul[1].header)

        # pre_source_table for collecting star info
        src_table = butler.get('preSourceTable',dataId={'instrument' : 'DECam', 'band':band, 'detector':det, 'visit':visit},collections='DECam/processing/calexp_{band}'.format(band=band))

        psf_candidates = src_table[src_table[psf_flag]]
    
        # compute psf fwhm
        visit_fwhm += psf_from_stars(calexp,psf_candidates)
    
        # collect detector ra/dec
        ra_dec_pairs.append( w.all_pix2world(np.array([[1024,2048]]),1,ra_dec_order=True)[0] )
    
        # I'm now done reading from calexp
        hdul.close()
    
        # load moments and compute ellipticities
        ixx = psf_candidates['ixx'].values
        iyy = psf_candidates['iyy'].values
        ixy = psf_candidates['ixy'].values
        
        e1 = (ixx - iyy)/(ixx + iyy)
        e2 = 2*ixy/(ixx + iyy)
        theta = np.arctan(2 * ixy / (ixx - iyy))/2
        
        # stokes parameters
        u = np.sqrt( e1**2 + e2**2 )*np.cos(theta)
        v = np.sqrt( e1**2 + e2**2 )*np.sin(theta)
        
        # now I have everything I need to update the visit summary
        visit_ra += psf_candidates['ra'].values.tolist()
        visit_dec += psf_candidates['dec'].values.tolist()
        
        visit_e1 += e1.tolist()
        visit_e2 += e2.tolist()
        
        visit_u += u.tolist()
        visit_v += v.tolist()
        
        visit_det += ((np.ones(len(psf_candidates)) * det ).tolist())
        psf_candidates_index += psf_candidates.index.tolist()
    
    e = time.time()
    
    # creating and saving the visit summary
    visit_summary = pd.DataFrame({'ra' : visit_ra, 'dec':visit_dec, 'fwhm':visit_fwhm, 'e1':visit_e1, 'e2':visit_e2, 'u':visit_u, 'v':visit_v, 'detector':visit_det},index=psf_candidates_index)
    
    visit_summary.to_csv("check_visit/summary_tables/{band}_visit_{exp}_summary.csv".format(exp=visit,band=band))
    
    # saving the ccd-coords
    ra_dec_pairs = np.array(ra_dec_pairs)
    ra_dec_df = pd.DataFrame(ra_dec_pairs,index=detectors,columns=['ra','dec'])
    ra_dec_df.to_csv("check_visit/summary_tables/{band}_visit_{exp}_wcs.csv".format(exp=visit,band=band))
    
    print("Process took {timer} seconds!".format(timer= int(e-s) ))
    return visit_summary
    
def psf_from_stars(calexp,psf_candidates,trim_width=10,plot=False):
    
    # array to record fwhm's
    fwhm_2d=[]
    
    # this fitting is taken from Shenming's original check_visit script
    #AE 09-2024: check for and mask NaN's first, then fit!
    for i in range(len(psf_candidates)):
        star_coord = [psf_candidates['x'].values[i], psf_candidates['y'].values[i]]
        star_coord = [int(np.round(coord)) for coord in star_coord]
        
        trim = calexp[star_coord[1]-trim_width:star_coord[1]+trim_width+1, star_coord[0]-trim_width:star_coord[0]+trim_width+1]
        mask = ~np.isnan(trim)
        
        # try/except to catch errors from fitting
        try:
            surf_init = models.Moffat2D(amplitude=trim[trim_width,trim_width], x_0=star_coord[0], y_0=star_coord[1], gamma=2., alpha=1.)
            fit_surf = fitting.LevMarLSQFitter()
        
            col_mat, row_mat = np.meshgrid(range(star_coord[0]-trim_width,star_coord[0]+trim_width+1),range(star_coord[1]-trim_width,star_coord[1]+trim_width+1))
            surf = fit_surf(surf_init, col_mat[mask], row_mat[mask], trim[mask])
        except:
            # This isn't a great fix by any means... but if it's just one bad star it shouldn't cause any issues :)
            print("Encountered an error during fitting, assinging 10px for one of the reference stars...")
            fwhm_2d.append(10)
            continue
            
        fwhm_2d.append(surf.fwhm)
        
    return fwhm_2d

if __name__ == '__main__':
    
    # checking for correct usage
    if len(sys.argv)!=3:
        print('Usage: python this_script.py band cores')
        sys.exit(1)
    
    band = sys.argv[1]
    cores = int(sys.argv[2])
    
    # create a butler to read/write data with
    butler = Butler('repo/repo',writeable=False)
    
    # collect the id of every dataset in calexp_band
    exposure_generator = butler.registry.queryDimensionRecords('exposure',datasets='calexp',collections='DECam/processing/calexp_{band}'.format(band=band))
    exposure_ids = []
    
    # iterate through the generator to record every exposure number
    for exposure in exposure_generator:
        exposure_ids.append(exposure.id)
    
    # each detector is a separate dataset, all of which have their records checked by the query
    # occassionally a detector will fail, so we query over all first detectors then select the unique ids
    # per-detector failures can be caused due to it being one of the bad detectors, holes in our catalogs, or corrupt data
    exposure_ids = np.unique(np.array(exposure_ids))
    
    print("Butler has found the following visits:")
    print(exposure_ids)
    
    # create a dictionary containing the good detectors per-exposure
    exposure_detectors = {}
    for id in exposure_ids:
        detectors = []
        
        # query for the detectors in a given visit
        detector_query = butler.registry.queryDimensionRecords('detector',datasets='calexp',collections='DECam/processing/calexp_{band}'.format(band=band),dataId={'instrument' : 'DECam' , 'visit' : id})
        
        # collect those detectors together and add them to the dictionary
        for result in detector_query:
            detectors.append(result.id)
        exposure_detectors[id] = detectors
        
    # I can only run up to 20 visits in parallel at a time...
    # so split the visits into (nearly) groups of 20
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
            arguments.append((visit,exposure_detectors[visit],'calib_psf_used',band))
            
        # create a pool for parallelization
        with Pool(len(process_visits)) as pool:
            output = pool.starmap(summary_from_stars, arguments)
    

