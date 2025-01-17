import numpy as np
from astropy.io import fits
import sys
import glob

# check for correct usage
# this catches BOTH any full crashes and exceptions, which is good bc either of those will cause ingestion to fail!

if len(sys.argv)!=4:
    print('Usage: python this_script.py check_this_file jobarray_dex jobarray_max')
    sys.exit(1)

array_dex = int(sys.argv[2])
array_num = int(sys.argv[3]) + 1 # start from zero go to N
files = glob.glob("raws/*.fits.fz")
chunk = np.array_split(files,array_num)[array_dex]

if sys.argv[1] not in chunk:
    print("That's not my problem!")
    quit()

check_me = fits.open(sys.argv[1])

# for some reason there aren't enough headers!
#if len(check_me) < 71:
#    raise Exception("Not enough headers, some CCDs are missing!")

# this opens each CCD to make sure there aren't errors in any particular header
# which could cause a segmentation fault!

for i in range(62):
    data = check_me[i+1].data

print("Data is intact!")

