import sys
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
import astropy.io.fits as fits
import astropy.utils as autils
from astropy.coordinates import SkyCoord
import requests
import json
from astropy.io import ascii
import subprocess
import glob
import time
import os

# this script starts and monitors the downloads for completion

# submit the first round of download_jobs based on the number of download_task_*.sh files that exist
# this let's me break the downloads into an arbitrary number of pieces for testing
completion_threshold=0
for i,job in enumerate(glob.glob("raws/download_task_*.sh")):
    completion_threshold+=1
    subprocess.run(['sbatch',f'raws/download_task_{i}.sh',str(i)])

# keep this while loop running until all the downloads finish!
while True:
    print("Loop starting!")
    counter = 0
    complete_flags = glob.glob("raws/*_download_complete*")
    print(complete_flags)
    print(len(complete_flags))
    # if a flag exists, check if all files were all downloaded
    if len(complete_flags) > 0:
        downloaded_files = glob.glob("raws/*.fits.fz")
        for flag in complete_flags:
        
            # lazy way of geting the download manager id which completed
            dl_num = flag.split('/')[-1]
            dl_num = dl_num.split('_')[0]
            
            # read the file to get the lines
            with open(f'raws/download_task_{dl_num}.sh','r') as f:
                lines = f.readlines()
            delete_lines = []
            
            # if this file is empty, update the counter and check the next file
            if len(lines) <= 7:
                print("file is empty!")
                counter+=1
                continue
            
            # make an array of line-numbers to delete
            for i, line in enumerate(lines):
                for filename in downloaded_files:
                    filename = filename.split('/')[-1]
                    if (filename in line) and ("curl" in line):
                        delete_lines.append(i-1)
                        delete_lines.append(i)
                        delete_lines.append(i+1)
                        break
            
            print(delete_lines)
            # now delete the lines with completed downloads
            with open(f'raws/download_task_{dl_num}.sh','w') as f:
            
                for i, line in enumerate(lines):
                    # skip the completed lines
                    if i in delete_lines:
                        print(f'skipping {line}')
                        continue
                    else:
                        print(f'writing {line}')
                        f.write(line)
            f.close()
            
            # delete the output file and resubmit the job
            print('deleting and resubmitting')
            os.remove(flag)
            subprocess.run(["sbatch",f'raws/download_task_{dl_num}.sh',str(dl_num)])
        
        # if all files have two or fewer lines, the downloads are finished!
        print('exited loop, checking counter')
        print(counter)
        if counter >= completion_threshold:
            break
        
    # wait 2min before checking the downloads again
    time.sleep(120) 
