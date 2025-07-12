import os
import time
import math
import glob
import sys

import numpy as np
import astropy as ap
import pandas as pd
import matplotlib.pyplot as pl

from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.visualization import ZScaleInterval
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

from astroquery.ipac.ned import Ned

from lsst.daf.butler import Butler

# create scatter-plot of fwhm against ellipticity
# not efficient to have two plotting functions here... but this only costs a few seconds per-band
def plot_scatter(visit_dict,fwhm_cut,ellip_cut,res=120):

    # statistics of interest
    fwhm_hist = []
    ellip_hist = []
    det_star = []

    # collect statistics per-detector
    visits = list(visit_dict)
    for visit in visits:
        for i in range(len(visit_dict[visit][1:])):
            det_star.append(visit_dict[visit][1:][i][3])
            fwhm_hist.append(visit_dict[visit][1:][i][5])
            ellip_hist.append(visit_dict[visit][1:][i][4])
    
    pl.scatter(ellip_hist,fwhm_hist,c=det_star,cmap='jet',alpha=0.5,marker='.')
    pl.axhline(fwhm_cut,color='Grey',linestyle='--')
    pl.axvline(ellip_cut,color='Grey',linestyle='--')
    pl.xlabel("Ellipticity")
    pl.ylabel("FWHM (px)")
    pl.title("FWHM v. Ellipticity")
    pl.colorbar()
    pl.savefig("check_visit/band_summary_plots/{band}_psf_ellip_scatter.png".format(band=band),dpi=res)
    
# create histogram across all detectors/visits
def plot_histograms(visit_dict,fwhm_cut,ellip_cut,good_visits=None,res=120):
    
    visits = list(visit_dict)
    
    # create histograms before cut
    if good_visits == None:
    
        # statistics of interest
        fwhm_hist = []
        ellip_hist = []
        det_star = []
        visit_star = []
    
        # collect statistics per-detector
        for visit in visits:
            visit_star.append(visit_dict[visit][0])
            for i in range(len(visit_dict[visit][1:])):
                det_star.append(visit_dict[visit][1:][i][3])
                fwhm_hist.append(visit_dict[visit][1:][i][5])
                ellip_hist.append(visit_dict[visit][1:][i][4])
        
        # create and save histograms+cuts
        pl.figure(1)
        pl.hist(det_star,bins='auto',histtype='step')
        pl.xlabel("Number of psf stars (calib_psf_used)")
        pl.ylabel("Count")
        pl.title("Distribution of PSF-Stars (per-detector)")
        pl.savefig("check_visit/band_summary_plots/{band}_psf_det_hist.png".format(band=band),dpi=res)
        pl.close(1)
        
        # create and save histograms+cuts
        pl.figure(2)
        pl.hist(fwhm_hist,bins='auto',histtype='step')
        pl.xlabel("FWHM (px)")
        pl.ylabel("Count")
        pl.title("Distribution of FWHM")
        pl.axvline(fwhm_cut,color='Grey',linestyle='--')
        pl.savefig("check_visit/band_summary_plots/{band}_fwhm_hist.png".format(band=band),dpi=res)
        pl.close(2)

        # create and save histograms+cuts
        pl.figure(3)
        pl.hist(ellip_hist,bins='auto',histtype='step')
        pl.xlabel("Ellipticity")
        pl.ylabel("Count")
        pl.title("Distribution of Ellipticities")
        pl.axvline(ellip_cut,color='Grey',linestyle='--')
        pl.savefig("check_visit/band_summary_plots/{band}_ellip_hist.png".format(band=band),dpi=res)
        pl.close(3)
        
        pl.figure(4)
        pl.hist(visit_star,bins='auto',histtype='step')
        pl.xlabel("Number of psf stars (calib_psf_used)")
        pl.ylabel("Count")
        pl.title("Distribution of PSF-Stars (per-visit)")
        pl.savefig("check_visit/band_summary_plots/{band}_psf_visit_hist.png".format(band=band),dpi=res)
        pl.close(4)
        
    else:
    
        # statistics of interest
        fwhm_hist = []
        ellip_hist = []
        det_star = []
        visit_star = []
    
        # collect statistics per-detector
        for visit in visits:
        
            # skip the bad visits
            if visit not in good_visits:
                continue
            
            visit_star.append(visit_dict[visit][0])
            for i in range(len(visit_dict[visit][1:])):
                
                # skip the bad detectors
                if visit_dict[visit][1:][i][0] not in good_visits[visit]:
                    continue
                    
                det_star.append(visit_dict[visit][1:][i][3])
                fwhm_hist.append(visit_dict[visit][1:][i][5])
                ellip_hist.append(visit_dict[visit][1:][i][4])
        
        # create and save histograms+cuts
        pl.figure(1)
        pl.hist(det_star,bins='auto',histtype='step')
        pl.xlabel("Number of psf stars (calib_psf_used)")
        pl.ylabel("Count")
        pl.title("Cut Distribution of PSF-Stars (per-detector)")
        pl.savefig("check_visit/band_summary_plots/{band}_cut_psf_det_hist.png".format(band=band),dpi=res)
        pl.close(1)
        
        # create and save histograms+cuts
        pl.figure(2)
        pl.hist(fwhm_hist,bins='auto',histtype='step')
        pl.xlabel("FWHM (px)")
        pl.ylabel("Count")
        pl.title("Cut Distribution of FWHM")
        pl.savefig("check_visit/band_summary_plots/{band}_cut_fwhm_hist.png".format(band=band),dpi=res)
        pl.close(2)

        # create and save histograms+cuts
        pl.figure(3)
        pl.hist(ellip_hist,bins='auto',histtype='step')
        pl.xlabel("Ellipticity")
        pl.ylabel("Count")
        pl.title("Cut Distribution of Ellipticities")
        pl.savefig("check_visit/band_summary_plots/{band}_cut_ellip_hist.png".format(band=band),dpi=res)
        pl.close(3)
        
        pl.figure(4)
        pl.hist(visit_star,bins='auto',histtype='step')
        pl.xlabel("Number of psf stars (calib_psf_used)")
        pl.ylabel("Count")
        pl.title("Cut Distribution of PSF-Stars (per-visit)")
        pl.savefig("check_visit/band_summary_plots/{band}_cut_psf_visit_hist.png".format(band=band),dpi=res)
        pl.close(4)
        
# this function actually applies the cuts we're interested in
# it returns visits in a dictionary followed by an array of the good detectors
def apply_cuts(visit_dict,visit_star_cut,det_star_cut,fwhm_cut,ellip_cut,det_cut=10,dist_cut=1.5):
    
    good_dict = {}
    for visit in list(visit_dict):
        print("Applying cuts to visit {exp}!".format(exp=visit))
        
        visit_statistics = visit_dict[visit]
        # if a visit does not pass the visit-level star-cut, throw it away
        if visit_statistics[0] < visit_star_cut:
            print("Visit {exp} failed the visit-level star-cut!".format(exp=visit))
            continue
        
        good_detectors = []
        # now apply detector-level cuts
        for i in range(len(visit_statistics[1:])):
            detector_statistics = visit_statistics[1:][i]
            
            # first a wcs cut, keep detectors w/in 1.5deg of the cluster center
            dist = cl_center.separation(SkyCoord(detector_statistics[1],detector_statistics[2],unit='deg')).deg
            if dist > 1.5:
                print("Visit {exp} det {num} failed the detector-level distance-cut!".format(exp=visit,num=detector_statistics[0]))
                continue
            
            # second a star-cut
            if detector_statistics[3] < det_star_cut:
                print("Visit {exp} det {num} failed the detector-level star-cut!".format(exp=visit,num=detector_statistics[0]))
                continue
            
            # third a fwhm cut
            if detector_statistics[5] > fwhm_cut:
                print("Visit {exp} det {num} failed the fwhm-cut!".format(exp=visit,num=detector_statistics[0]))
                continue
            
            # finally an ellipticity cut
            if detector_statistics[4] > ellip_cut:
                print("Visit {exp} det {num} failed the ellip-cut!".format(exp=visit,num=detector_statistics[0]))
                continue
            
            # detectors reaching here match our quality-cuts
            good_detectors.append(detector_statistics[0])
        
        # if there are very few good detectors, we throw away the exposure
        if len(good_detectors) <= det_cut:
            print("Visit {exp} only has {num} < {cut} good exposures, throwing it away!\n".format(exp=visit,num=len(good_detectors),cut=det_cut))
            continue
        else:
            print("Visit {exp} has {num} good exposures, so we'll keep them. The good detectors are:".format(exp=visit,num=len(good_detectors)))
            print(good_detectors)
            print()
            good_dict[visit] = good_detectors
    
    return good_dict

if __name__ == '__main__':
    
    # checking for correct usage
    if len(sys.argv)!=5:
        print('Usage: python this_script.py cln band fwhm_px ellip')
        sys.exit(1)
    
    cluster_name = sys.argv[1]
    band = sys.argv[2]
    fwhm_cut = float(sys.argv[3])
    ellip_cut = float(sys.argv[4])
    
    # load info about the cluster from NED
    ned_result = Ned.query_object(cluster_name)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    
    cl_center = SkyCoord(ra_cl, dec_cl, unit="deg")
    
    # create a butler to read/write data with
    butler = Butler('repo/repo',writeable=True)
    
    
    # load visits from the butler (reused from check_visit, not optimal but querying is fast
    exposure_generator = butler.registry.queryDimensionRecords('exposure',datasets='calexp',collections='DECam/processing/calexp_{band}'.format(band=band))
    exposure_ids = []
    
    for exposure in exposure_generator:
        exposure_ids.append(exposure.id)
    
    exposure_ids = np.unique(np.array(exposure_ids))
    
    # collect the visit/detector statistics in a dictionary at the visit-level
    # each entry in the dict is an array, first entry is total # of psf-stars, following entries are per-det arrays which store:
    # [detector #, ra_center, dec_center, num-psf-stars, mean-ellip, mean-fwhm]
    # I'll also create separate arrays for counting per detector/visit psf_stars
    # there is probably a better way of doing this... but this works and its fast :)
    
    visit_dict = {}
    per_visit_psf_stars = []
    per_det_psf_stars = []
    
    for visit in exposure_ids:
        fwhm_cut
        visit_statistics = []
        
        # loading the visit summary
        try:
            visit_summary = pd.read_csv("check_visit/summary_tables/{band}_visit_{exp}_summary.csv".format(band=band,exp=visit),index_col=0)
        except:
            print(f'Cant find visit summary for band {band} exposure {visit}, skipping it!')
            continue
            
        # loading the det wcs
        visit_wcs = pd.read_csv("check_visit/summary_tables/{band}_visit_{exp}_wcs.csv".format(band=band,exp=visit),index_col=0)
        
        # now we'll filter-out the bad detectors
        visit_summary_bad = (visit_summary['detector'] == 2) | (visit_summary['detector'] == 31) | (visit_summary['detector'] == 61) | (visit_summary['detector'] == 62)
        visit_wcs_bad = (visit_wcs.index == 2) | (visit_wcs.index == 31) | (visit_wcs.index == 61) | (visit_wcs.index == 62)
    
        visit_summary = visit_summary.drop(visit_summary[visit_summary_bad].index)
        visit_wcs = visit_wcs[np.invert(visit_wcs_bad)]
        
        # append the total number of psf_used stars
        visit_statistics.append(len(visit_summary))
        
        # append this to a separate array for computing visit mean/std
        per_visit_psf_stars.append(len(visit_summary))
        
        # collect detector statistics
        detectors = np.unique(visit_summary['detector'])
        for det in detectors:
            
            detector_summary = visit_summary[visit_summary['detector'] == det]
            det_statistics = []
            det_statistics.append(int(det))
            
            det_statistics.append(visit_wcs['ra'][det])
            det_statistics.append(visit_wcs['dec'][det])
    
            det_statistics.append(len(detector_summary))
            
            # appending this to a separate array for computing det mean/std
            per_det_psf_stars.append(len(detector_summary))
            
            det_statistics.append(np.mean( np.sqrt(detector_summary['e1']**2 + detector_summary['e2']**2)))
            det_statistics.append(np.mean( detector_summary['fwhm']))
            
            visit_statistics.append(det_statistics)
        
        visit_dict[visit] = visit_statistics
    
    visit_star_cut = np.mean(per_visit_psf_stars) - 2*np.std(per_visit_psf_stars)
    det_star_cut = np.mean(per_det_psf_stars) - 3*np.std(per_det_psf_stars)
    
    # I could parallelize this function... but it really doesn't take that long to begin with so I don't think its worthwhile
    good_visits = apply_cuts(visit_dict,visit_star_cut,det_star_cut,fwhm_cut=fwhm_cut,ellip_cut=ellip_cut)
    
    # create histogram before and after quality-cuts 
    print("Plotting pre and post-cut histograms")
    plot_histograms(visit_dict,fwhm_cut,ellip_cut,good_visits=None,res=120)
    plot_histograms(visit_dict,fwhm_cut,ellip_cut,good_visits,res=120)
    plot_scatter(visit_dict,fwhm_cut,ellip_cut,res=120)
    
    # now its time to write these to the butler
    print("Creating new collection containing the good detectors: DECam/processing/quality_detectors_{band}".format(band=band))
    
    # these are the dataset types required later in the pipeline
    # might need to add more depending on what later steps require...
    dataset_types=['src','calexp','preSourceTable','src_schema','calexpBackground']
    
    # formatting the query, sloppy but it works
    where_string = "instrument='DECam' AND band='{filt}' AND ( ".format(filt=band)
    
    for visit in list(good_visits):
        
        where_string += " ( (visit = {exp}) AND (".format(exp=visit)
        for det in good_visits[visit]:
            where_string += " (detector={det}) OR".format(det=det)
        where_string = where_string[:-2] + ") ) OR"
    
    where_string = where_string[:-2] + ")"
    
    # register the collection and associate the datasets
    butler.registry.registerCollection("DECam/processing/quality_detectors_{band}".format(band=band))
    refs = butler.registry.queryDatasets(dataset_types,collections='DECam/processing/calexp_{band}'.format(band=band),where=where_string,findFirst=True)
    butler.registry.associate("DECam/processing/quality_detectors_{band}".format(band=band),refs)
    
    # eventually DM will move way from registry, when/if it becomes depricated the code-snippet below will collect the datasets
    #for dtype in dataset_types:
        #refs = butler.query_datasets(dtype,collections='DECam/processing/calexp_{band}'.format(band=band),where=where_string,find_first=True)
        #butler.registry.associate("DECam/processing/quality_detectors_{band}".format(band=band),refs)

    
