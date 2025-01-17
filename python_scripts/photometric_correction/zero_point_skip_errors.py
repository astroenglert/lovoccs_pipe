# zero-point correction
# apply zero-point corrections saved in the residual_table

import sys
import os

from astropy.table import Table

import numpy as np

#TODO this should really be merged with zero_point.py
#homebrew modules here, these will likely change in the future
from color_terms import get_instrument_headers

if __name__ == '__main__':
    
    if len(sys.argv) != 5:
        print("python zero_point.py catalog color_term_residual stellar_locus_residual output_filename")
        raise Exception("Improper Usage! Correct usage: python zero_point.py catalog color_term_residual stellar_locus_residual output_filename")
    
    # loading arguments from cln
    catalog_filename = sys.argv[1]
    ct_bias_filename = sys.argv[2]
    sl_bias_filename = sys.argv[3]
    output_filename = sys.argv[4]
    
    # loading tables from disk
    catalog = Table.read(catalog_filename,format='ascii.csv')
    ct_bias = Table.read(ct_bias_filename,format='ascii.csv')
    sl_bias = Table.read(sl_bias_filename,format='ascii.csv')
    
    # loading magnitudes which exist for decam
    # to apply this to a different instrument, this will have to be adjusted
    instr = 'decam'
    headers = get_instrument_headers(instr)
    column_names = catalog.colnames
    magnitude_list = []
    for entry in list(headers.values()):
        if (entry[-4:] == '_mag') & (entry in column_names):
            magnitude_list.append(entry)
    
    # now start applying the corrections!
    err_tag = 'err'
    for band in magnitude_list:
        print(f'correcting {band}')
        # use ct_correction for everything EXCEPT u-band and Y-band (TODO manually skip the latter for now... )
        if band[0] == 'Y':
            catalog[band] = catalog[band]
        elif band[0] == 'u':
            catalog[band] = catalog[band] - sl_bias[band[0]][0]
            #catalog[band + err_tag] = np.sqrt(catalog[band + err_tag]**2 + sl_bias[band[0]][1]**2)
        else:
            catalog[band] = catalog[band] - ct_bias[band[0]][0]
            #catalog[band + err_tag] = np.sqrt(catalog[band + err_tag]**2 + ct_bias[band[0]][1]**2)
        
    catalog.write(output_filename, format="ascii.csv", overwrite=True)
    
    
    
    
    
    


