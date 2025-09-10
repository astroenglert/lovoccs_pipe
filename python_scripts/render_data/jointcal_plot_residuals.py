import glob
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from astropy.io import ascii
from astropy.coordinates import search_around_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

def photom_astrom_jc_match(photom_table,astrom_table):
    '''
    Helper function to load the results of jointcal debug tables into a more managable format
    
    Args:
      photom_table: DataFrame; dataframe containing jointcal photometry output
      astrom_table: DataFrame; dataframe containing jointcal astrometry output
    
    Returns:
      output: Astropy Table; a table containing the essential columns for debugging jointcal
    
    '''
    photom_ids = photom_table['id'].astype(int)
    astrom_ids = astrom_table['id'].astype(int)
    shared_astrom = np.in1d(astrom_ids,photom_ids)
    shared_photom = np.in1d(photom_ids,astrom_ids)
    photom_table = photom_table[shared_photom]
    astrom_table_shared = astrom_table[shared_astrom]
    star_id = astrom_table_shared['id'].astype(int)
    ra = astrom_table_shared['ra'].astype(float)
    dec = astrom_table_shared['dec'].astype(float)
    fitted_mag = photom_table['mag'].astype(float)
    inst_mag = photom_table['instMag'].astype(float)
    inst_flux = photom_table['instFlux'].astype(float)
    input_flux = photom_table['inputFlux'].astype(float)
    fitted_flux = photom_table['fittedFlux'].astype(float)
    num_meas = photom_table['nm'].astype(float)
    det_id = photom_table['detector'].astype(int)
    visit_id = photom_table['visit'].astype(int)
    output =  Table(data={"ra" : ra,"dec" : dec,"fittedMag" : fitted_mag,"instMag" : inst_mag,"insFlux": inst_flux,"inputFlux":input_flux,"fittedFlux":fitted_flux,"nm":num_meas,"detector":det_id,"visit":visit_id,"star_id":star_id})
    return output

# load a reference catalog
def get_reference_catalog(refcat_name,cluster_name,band=None):
    '''
    Helper function to load a reference catalog
    
    Args:
      refcat_name: string; name of refcat to load
      cluster_name: string; name of the cluster
    
    Returns:
      return_me: Table; astropy table storing reference stars
    
    '''
    
    # load the reference catalog
    cluster_refcats = glob.glob("/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new/{CLN}/{REFCAT}*.csv".format(CLN=cluster_name,REFCAT=refcat_name))
    
    if len(cluster_refcats) == 0:
        print("Uh-oh, I couldn't find any catalogs matching that name!")
        raise Exception()
    
    # sm is a little special so we need to handle it differently
    if len(cluster_refcats) > 1 and "sm" not in refcat_name:
        print("Uh-oh, you weren't specific enough in your choice of reference catalog, I received multiple potential catalogs matching that description")
        print(cluster_refcats)
        raise Exception()
    
    ra_header = "ra"
    dec_header = "dec"
    
    # filter-maps and coordinate headers
    if "des" in refcat_name:
        mag_col = "wavg_mag_psf_" + band
    
    if "legacy" in refcat_name:
        mag_col = "mag_" + band
    
    if "ps1" in refcat_name:
        mag_col = band + "mag"
        ra_header = "RAJ2000"
        dec_header = "DEJ2000"
        
    if "sdss" in refcat_name:
        if band != "u":
            print("Uh-oh, we only have u-band for SDSS!")
            raise Exception()
        mag_col = "upmag"
        ra_header = "RA_ICRS"
        dec_header = "DE_ICRS"
        
    if "sm" in refcat_name:
        if band == "u":
            band_search = "v"
            mag_col = "v_psf"
        else:
            band_search = band
            mag_col = band + "_psf"
        ra_header = "raj2000"
        dec_header = "dej2000"
        
        # sm has different filename when we download them, so its a little special
        cluster_refcats = glob.glob("/gpfs/data/idellant/Clusters/calib_catalog_repo/catalogs_new/{CLN}/{REFCAT}*_{CLN}_{BAND}.csv".format(CLN=cluster_name,REFCAT=refcat_name,BAND=band_search)) 
        
        # this should return a single catalog, if not something is wrong!
        if len(cluster_refcats) != 1:
            print("Uh-oh, I found more sm catalogs than I should've!")
            print(cluster_refcats)
    
    # for jointcal, I only really care about magnitude and ra/dec
    table = ascii.read(cluster_refcats[0],header_start=0,delimiter=",")
    return_me = table[ra_header,dec_header,mag_col]
    return_me.rename_column(mag_col,band)
    return_me.rename_column(ra_header,'ra')
    return_me.rename_column(dec_header,'dec')
    
    return return_me

def compare_refcat_jointcal(cluster_name,refcat_name,band,order,jointcal_dir=None,visit=None,magu=25,magl=0,cmin=-0.2,cmax=0.2):
    '''
    Plot the jointcal residuals
    
    Args:
      cluster_name: string; name of the cluster
      refcat_name: string; name of the refcat
      band: string; band to load
      order: int; order of the jointcal fit
      jointcal_dir: string; string to save residual plots in, defaults to .../jointcal_residuals/
      visit: int; visit-id to draw
      magu: float; upper magnitude limit to check
      magl: float; lower magnitude limit to check
      cmin: float; minimum residual on colorbar
      cmax: float; maximum residual on colorbar
    
    Returns:
      None
    
    '''
    
    # will need to standardize these filepaths
    final="final"
    
    # by default it will search jointcal_residuals for residuals of the given-order
    if jointcal_dir == None:
        jointcal_photom_path = glob.glob("jointcal_residuals/order_{ORD}/photometry_{FINAL}*_{BAND} *meas.csv".format(ORD=order,FINAL=final,BAND=band))
        jointcal_astrom_path = glob.glob("jointcal_residuals/order_{ORD}/astrometry_{FINAL}*_{BAND} *meas.csv".format(ORD=order,FINAL=final,BAND=band))
    else:
    # but if the files are stored somewhere else it will search for them in that directory

        jointcal_photom_path = glob.glob(jointcal_dir + "photometry_{FINAL}*_{BAND} *meas.csv".format(FINAL=final,BAND=band))
        jointcal_astrom_path = glob.glob(jointcal_dir + "astrometry_{FINAL}*_{BAND} *meas.csv".format(FINAL=final,BAND=band))
    
    if len(jointcal_photom_path) != 1 or len(jointcal_astrom_path) != 1:
        print("Uhoh, I can't find the jointcal debug catalogs")
        raise Exception()
    
    jointcal_photom_path = jointcal_photom_path[0]
    jointcal_astrom_path = jointcal_astrom_path[0]

    # read all of the data into python
    ext_survey = get_reference_catalog(refcat_name,cluster_name,band=band)
    photom_table = ascii.read(jointcal_photom_path,delimiter='\t')
    astrom_table = ascii.read(jointcal_astrom_path,delimiter='\t')
    photom_table = photom_table[1:] #first rows are comments
    astrom_table = astrom_table[1:] 
    
    # apply magnitude cuts
    photom_table = photom_table[(photom_table['mag'].astype(float).value < magu) & (photom_table['mag'].astype(float).value > magl)]

    # pick-out a specific visit if one is specified
    if visit != None:
        photom_table = photom_table[photom_table['visit'].astype(int) == visit]
        astrom_table = astrom_table[astrom_table['visit'].astype(int) == visit]

    # use LSP star-id to collect jc-coords+mags
    matched = photom_astrom_jc_match(photom_table,astrom_table)

    # remove the reference sources
    matched = matched[matched['nm'] > 1]

    # load coordinates from external catalog+jc and run source matching
    ext_coords = SkyCoord(ra=ext_survey['ra']*u.deg,dec=ext_survey['dec']*u.deg)
    jc_coords = SkyCoord(ra=matched['ra']*u.deg,dec=matched['dec']*u.deg)

    idjc,idext,d2d,d3d = search_around_sky(jc_coords,ext_coords,1*u.arcsec)
    jc_matches = matched[idjc]
    ext_matches = ext_survey[idext]

    # all rows are matched at this point, so plotting from here is easy!

    # first, lets create a scatter-plot for the magnitudes relative to the reference catalog
    pl.figure()

    pl.scatter(jc_matches['ra'],jc_matches['dec'],c=jc_matches['fittedMag'] - ext_matches[band],cmap='jet',alpha=0.4,marker='.',vmin=cmin,vmax=cmax,s=3)
    pl.colorbar()
    pl.xlabel('RA (deg)')
    pl.ylabel('DEC (deg)')

    pl.title("JC-EXT Residuals")
    
    # these plots are internal for now, we can fix them up if we want to later on
    if visit == None:
        pl.savefig("jointcal_residuals/{BAND}_jointcal-{REFCAT}_o{ORD}_residual.png".format(ORD=order,REFCAT=refcat_name,BAND=band))
    else:
        pl.savefig("jointcal_residuals/{BAND}_jointcal-{REFCAT}_o{ORD}_residual_visit_{VISIT}.png".format(ORD=order,REFCAT=refcat_name,VISIT=visit,BAND=band))

    pl.close()
    
    # next lets make a histogram of the residuals (from before and after lsq-opt), first I need to collect the initial photom
    
    final="initial"
    
    # by default it will search jointcal_residuals of A85 for residuals of the given-order (just from testing)
    if jointcal_dir == None:
        jointcal_photom_path = glob.glob("/gpfs/data/idellant/englert_newPipelineDev/A85/jointcal_residuals/order_{ORD}/photometry_{FINAL}*_{BAND} *meas.csv".format(ORD=order,FINAL=final,BAND=band))
        jointcal_astrom_path = glob.glob("/gpfs/data/idellant/englert_newPipelineDev/A85/jointcal_residuals/order_{ORD}/astrometry_{FINAL}*_{BAND} *meas.csv".format(ORD=order,FINAL=final,BAND=band))
    else:
    # if the files are stored somewhere else it will search for them in that directory
        jointcal_photom_path = glob.glob(jointcal_dir + "photometry_{FINAL}*_{BAND} *meas.csv".format(FINAL=final,BAND=band))
        jointcal_astrom_path = glob.glob(jointcal_dir + "astrometry_{FINAL}*_{BAND} *meas.csv".format(FINAL=final,BAND=band))
    
    if len(jointcal_photom_path) != 1 or len(jointcal_astrom_path) != 1:
        print("Uhoh, I can't find the jointcal debug catalogs")
        raise Exception()
    
    jointcal_photom_path = jointcal_photom_path[0]
    jointcal_astrom_path = jointcal_astrom_path[0]

    # read all of the data into python
    ext_survey = get_reference_catalog(refcat_name,cluster_name,band=band)
    photom_table = ascii.read(jointcal_photom_path,delimiter='\t')
    astrom_table = ascii.read(jointcal_astrom_path,delimiter='\t')
    photom_table = photom_table[1:] #first rows are comments
    astrom_table = astrom_table[1:] 
    
    # apply magnitude cuts
    photom_table = photom_table[(photom_table['mag'].astype(float).value < magu) & (photom_table['mag'].astype(float).value > magl)]

    # pick-out a specific visit if one is specified
    if visit != None:
        photom_table = photom_table[photom_table['visit'].astype(int) == visit]
        astrom_table = astrom_table[astrom_table['visit'].astype(int) == visit]

    # use LSP star-id to collect jc-coords+mags
    matched = photom_astrom_jc_match(photom_table,astrom_table)

    # remove the reference sources
    matched = matched[matched['nm']>1]

    jc_coords = SkyCoord(ra=matched['ra']*u.deg,dec=matched['dec']*u.deg)

    idjc,idext,d2d,d3d = search_around_sky(jc_coords,ext_coords,1*u.arcsec)
    jc_matches_initial = matched[idjc]
    ext_matches_initial = ext_survey[idext]
    
    final_residuals = jc_matches['fittedMag'] - ext_matches[band]
    initial_residuals = jc_matches_initial['fittedMag'] - ext_matches_initial[band]
    
    fig,ax = pl.subplots()
    pl.hist(final_residuals,bins=100,label="post-jointcal",histtype='step')
    pl.hist(initial_residuals,bins=100,label="pre-jointcal",histtype='step')
    pl.legend()
    
    # need to settle on scaling
    pl.xlim((-0.5,0.5)) # this encloses almost all of the residuals except outliers
    pl.xlabel("$ m_{jc} - m_{ref} $")
    pl.ylabel("Count")
    pl.title("Residual Histogram")
    
    # I need to compute the mean and, more importantly, the stdv of this distribution
    
    mean = np.mean(final_residuals)
    stdv = np.std(final_residuals)
    
    mean_initial = np.mean(initial_residuals)
    stdv_initial = np.std(initial_residuals)
    
    pl.text(x=0.2,y=0.9,s="$ \\mu = {0:.3f} $ \n $ \\sigma = {1:.3f} $".format(mean,stdv),transform=ax.transAxes,ha='center')
    
    pl.axvline(mean,c='C0',ls='--')
    pl.axvline(mean_initial,c='C1',ls='--')
    
    if visit == None:
        pl.savefig("jointcal_residuals/{BAND}_jointcal-{REFCAT}_o{ORD}_residual_hist.png".format(ORD=order,REFCAT=refcat_name,BAND=band))
    else:
        pl.savefig("jointcal_residuals/{BAND}_jointcal-{REFCAT}_o{ORD}_residual_visit_{VISIT}_hist.png".format(ORD=order,REFCAT=refcat_name,VISIT=visit,BAND=band))
    
    pl.close()
    
    pl.figure()
    
    pl.plot(ext_matches[band],final_residuals,'.',alpha=0.4,ms=3)
    pl.xlabel("$ m_{ref} $")
    pl.ylabel("$ m_{jc} - m_{ref} $")
    pl.title("Residual v. $ m_{ref} $")
    
    if visit == None:
        pl.savefig("jointcal_residuals/{BAND}_jointcal-{REFCAT}_o{ORD}_residual_plot.png".format(ORD=order,REFCAT=refcat_name,BAND=band))
    else:
        pl.savefig("jointcal_residuals/{BAND}_jointcal-{REFCAT}_o{ORD}_residual_visit_{VISIT}_plot.png".format(ORD=order,REFCAT=refcat_name,VISIT=visit,BAND=band))

    pl.close()

# update this to run after step 2b

# these functions are really useful to have as a debug tool in python...
# this ensures that if the script is imported the code below here will not run
if __name__ == '__main__':
    
    if len(sys.argv)!=5:
        raise Exception('Usage: python jointcal_residuals.py cln jointcal_dir refcat order band')

    cluster_name = sys.argv[1]
    refcat = sys.argv[2]
    order = sys.argv[3]
    band = sys.argv[4]
    
    if refcat == "sdss":
        refcat = "sdss_dr12_{CLN}".format(CLN=cluster_name)
    
    # by default these files are saved in {CLN}/jointcal_residuals
    compare_refcat_jointcal(cluster_name,refcat,band,jointcal_dir="jointcal_residuals/order_{ORD}/".format(ORD=order),order=order,cmin=-0.5,cmax=0.5)
    
    
    

