from astropy.table import Table
import numpy as np
import os
import sys
from . import utilities
from .utilities import new_columns

#TODO is it a good practice to override catalog values here to zero? How does that impact/change what the calibration does? I thought it used statistics of the entire sample and if we're padding it with zeros that may not be good.....
def fix_nan(catalog, key):
    """
    Routine to fix NaN entries.
    """
    x = catalog[key]
    mask = np.isfinite(x)
    n_fix = mask.astype(int).sum()
    if n_fix > 0:
        catalog[key][mask] = 0.

def main(argv):
    """
    Main function to output a FITS catalog with calibrations and weights, given some input catalog
    based on the LSST Science Pipelines image processing routines.

    Usage: The routine is run as follows:

        python gen_hsc_calibrations.py input_file output_file_name

    It will read in the input file (`input_file`) in FITS format, and produce an output FITS file
    named `output_file_name` with the new information.

    Dependencies: To run this script, you will need to have installed astropy and numpy.

    Memory usage: The script will ingest the data into memory, and make some additional columns, so
    the memory usage can range from 1-2x the size of the ingested catalog.
    """
    # Sanity check the arguments.
    if len(argv) != 3:
        raise RuntimeError("Wrong number of args; see below for usage info!\n"+main.__doc__)

    input_file = sys.argv[1]
    output_file_name = sys.argv[2]

#    if os.path.exists(output_file_name):
#        raise IOError("Output file %s already exists!  Exiting..."%output_file_name)

    if not os.path.exists(input_file):
        raise IOError("Input file %s does not exist!  Exiting..."%input_file)
    else:
        catalog = Table.read(input_file)
        print("Read in catalog with %d entries from %s"%(len(catalog),input_file))

    # Basic sanity checks: does the catalog have the entries needed?
    try:
        utilities.get_snr(catalog)
    except:
        raise RuntimeError("Input catalog does not have information needed to get cmodel SNR!")
    try:
        utilities.get_res(catalog)
    except:
        raise RuntimeError("Input catalog does not have information needed to get resolution factor!")
    try:
        utilities.get_psf_ellip(catalog)
    except:
        raise RuntimeError("Input catalog does not have PSF shape information, needed for additive shear bias corrections!")
        
    # Create the new entries.
    print("Creating weights, calibrations, etc.")
    calibration_info = utilities.create_calibs(catalog)

    # Fix NaN/inf values, which could result from some NaN/inf values in the input catalog.
    # Typically there are very few of these.
    for entry in new_columns:
        fix_nan(calibration_info, entry)

    print('Writing to file %s'%output_file_name)
    calibration_info.write(output_file_name, overwrite=True)

if __name__ == "__main__":
    main(sys.argv)
