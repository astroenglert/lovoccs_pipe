# Get raws and mastercals from noao(noirlab)
# https://github.com/NOAO/nat-nb
#print('Usage: python this_script.py cluster_name download_command_filename')


#=======================
import requests
import json
import sys
import random
from time import sleep

import numpy as np

import pandas as pd

import matplotlib.pyplot as pl

import astropy.units as u
from astropy.coordinates import get_body
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.coordinates import ICRS
from astropy.coordinates import EarthLocation
from astropy.coordinates import uniform_spherical_random_surface
from astropy.time import Time

from astroquery.ned import Ned

# I should have been using connection_pooling/keepAlive from the very beginning!
# this keeps a connection to astroarchive.noirlab.edu open... so we don't have to reconnect
# and get turned-down with a 500 response code!
session = requests.Session()

#=======================
def prompt(string):
    print('\n'+'-'*25+'\n'+string+'\n'+'-'*25)

# filter map for qol
filter_map = {
            'u': 'u DECam c0006 3500.0 1000.0',
            'g': 'g DECam SDSS c0001 4720.0 1520.0',
            'r': 'r DECam SDSS c0002 6415.0 1480.0',
            'i': 'i DECam SDSS c0003 7835.0 1470.0',
            'z': 'z DECam SDSS c0004 9260.0 1520.0',
            'Y': 'Y DECam c0005 10095.0 1130.0',
            'VR': 'VR DECam c0007 6300.0 2600.0',
    }

# generic function to query noirlab for raws
# derived from my query built for flats originally
def query_raws(caldat_min, caldat_max, ra, dec,
                # outfields, headers of output table
                outfields=["ra_center","dec_center","exposure","dateobs_center","caldat","ifilter","archive_filename","url","proposal"],
                # searchfields, format here is [ ["category1","element1"] , ["category2","element2"] , ... ]
                #TODO should I worry about the release date?
                search= [ ["instrument","decam"],["obs_type","object"],["proc_type","raw"],["prod_type","image"], ],
                # by default, do not limit the search results
                limit=None,
                ifilter=None,
                # by default we search in a 4degx4deg box, encloses (almost) all pointings w/ cluster
                radius=2,
                expmin=10,
                expmax=1000,
                session=session,
               ):
    '''
    Function to query for raw exposures from NOIRLab
    
    Args:
      caldat_min: string; the minimum calendar date to query
      caldat_maxd: string; the maximum calendar date to query
      ra: float; the Ra (in degrees) to center the query on
      dec: float; the Dec (in degrees) to center the query on
      outfields: array; list of fileds the query should return, 
                        defaults to ["ra_center","dec_center","exposure","dateobs_center","caldat","ifilter","archive_filename","url","proposal"]
      search: array; list specifying criteria for specific columns,
                     defaults to [ ["instrument","decam"],["obs_type","object"],["proc_type","raw"],["prod_type","image"], ]
      limit: int; maximum number of search results to return
      ifilter: string; the filter name to query
      radius: float; radius in deg to query about
      expmin: float; minimum exposure time
      expmax: float; maximum exposure time
      session: Session; Session for the connection
    
    Returns:
      dfMaster: DataFrame; dataframe containing visit information queried
    
    '''
    
    # roots for querying the server
    search.append(["caldat",caldat_min,caldat_max])
    natroot="https://astroarchive.noirlab.edu"
    adsurl = '%s/api/adv_search'%natroot
    apiurl = f'{adsurl}/find/?limit={limit}'
    
    # search within a radius-deg box of the specified pointing
    search.append(['ra_center',ra - radius/np.cos(dec*np.pi/180),ra + radius/np.cos(dec*np.pi/180)])
    search.append(['dec_center',dec - radius,dec + radius])
    
    # appending constraints on exptime we use for LV (60s-180s)
    search.append(['exposure',expmin,expmax])
    
    # if a filter is specified, append it to the search
    if ifilter is not None:
        search.append(["ifilter",filter_map[ifilter]])
    
    # default format for queries is a json, so lets build one from the input search
    json = { "outfields" : outfields , "search" :  search }
    
    print("Submitting query to NoirLab...",flush=True)
    
    # occassionally a json isn't returned from the query... server-side error likely so force a loop until that error doesn't occur!
    while True:
        try:
            request = session.post(apiurl,json=json)
            if request.reason == "Too Many Requests":
                print("Uh-oh we've sent too many requests!",flush=True)
                seconds = int(request.headers['Retry-After'])

                # wait for requested time then try again
                print("Lets wait " + str(seconds) + " seconds before continuing...",flush=True)
                sleep( seconds + random.random() * 0.1 )

                # resuming sending requests after the wait
                request = session.post(apiurl,json=json)

                # checking again for bad gateway then trying
                while request.reason == "Bad Gateway":
                    print("Bad Gateway! Trying again in 5s.",flush=True)
            
                    sleep( 5 + random.random() * 0.1 )
            
                    # this shouldn't cause a 429 Error due to the delay
                request = session.post(apiurl,json=json)

            dfMaster = pd.DataFrame(request.json()[1:])
        except:
            print("Encountered an error... here is the output request from Noirlab",flush=True)
            print(request,flush=True)
            print(request.content,flush=True)
            sleep( 15 + random.random() * 5 )
            continue
        break

    return dfMaster

def make_download_config(table=None,urls=None,filename=None,download_config='raws/download_task'):
    '''
    Write a file to execute a list of download commands
    
    Args:
      table: DataFrame; table containing urls and filenames to download
      urls: array; a list of URLs
      filename: array; a list of filenames corresponding to each URL
      download_config: string; tag to pre-pend to each script
    
    Returns:
      None 
    
    '''
    # collect url and archive filenames
    if table is not None:
        print("Pulling urls and filenames from the table.")
        urls = table['url']
        filename = table['archive_filename']
        iter_me = table.index
        
    else:
        iter_me = range(len(urls))
    
    print(iter_me)
    # download config file
    # create 8 files, store them in an array and cycle through them appending commands
    dl_files = []
    for i in range(8):
        dl_files.append(open(download_config + f'_{i}.sh','a'))
        
        # delay the start of each script just a little and add preamble for slurm
        dl_files[i].write('#!/bin/sh\n')
        dl_files[i].write('#SBATCH --time=24:00:00\n')
        dl_files[i].write(f'#SBATCH -o slurm_outputs/download_batch_{i}-%j.err\n')
        dl_files[i].write(f'#SBATCH -o slurm_outputs/download_batch_{i}-%j.out\n')
        dl_files[i].write('cd raws\n')
        dl_files[i].write(f'sleep {i}s\n')
        
    for i in range(len(iter_me)):
        index = iter_me[i]
        
        # need stricter restart settings, 100k/s is a sign that the connection is throttling
        # so restart with 5s below that after a 30s delay (to hopefully get to another port)
        # retry downloads 5 times with 30s delay
        # and wait 5s between downloads to avoid conflicts/resource limits
        
        dl_files[i%8].write('echo Now downloading %s\n'%(filename[index].split('/')[-1]))
        line = 'curl --remove-on-error --retry 0 --retry-delay 5 --retry-connrefused --speed-limit 200000 --speed-time 5 -o %s %s\n'%(filename[index].split('/')[-1], urls[index])
        dl_files[i%8].write(line)
        dl_files[i%8].write('sleep 5s\n')
        
    # close the txt file
    for i,file in enumerate(dl_files):
        file.write(f'echo I finished downloading! > {i}_download_complete.txt')
        file.close()

if __name__ == '__main__':
    
    #=======================
    # Pre-process
    if len(sys.argv)!=2:
        print('Usage: python this_script.py cluster_name')
        sys.exit(1)
    
    cluster_name = sys.argv[1]
    
    prompt('NED query...')
    ned_result = Ned.query_object(cluster_name)
    
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    redshift_cl = ned_result[0]['Redshift']
    print('%s: ra_cl, dec_cl, redshift_cl: '%cluster_name, ra_cl, dec_cl, redshift_cl)
    
    #TODO update caldate-max once calibrations are done
    df = query_raws('2012-01-01','2025-07-30',ra_cl,dec_cl)
    df.to_csv('query_result_%s.csv'%cluster_name)
    
    make_download_config(table=df)

