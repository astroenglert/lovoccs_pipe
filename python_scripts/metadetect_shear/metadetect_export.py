
import os
import time
import math
import glob
import sys
import json

from importlib import resources as impresources
from pathlib import Path

import numpy as np

from astropy.table import Table, vstack

from lsst.daf.butler import Butler

# homebrew modules below
from ..export.export_data import export_patch_data
from . import export_config
export_config = impresources.files(export_config)

if __name__ == '__main__':
    
    # checking for correct usage
    if len(sys.argv)!=2:
        print("python this_script.py cln")
        sys.exit(1)
    
    cln = sys.argv[1]
    
    # create a butler to read/write data with
    butler = Butler('repo/repo')
    
    # load the flags and columns
    flag_flp = export_config.joinpath('object_flags.json')
    cols_flp = export_config.joinpath('export_columns.json')
    
    with open(flag_flp) as f:
        flags = json.load(f)
    with open(cols_flp) as f:
        columns = json.load(f)
    
    # collect the patch indices for iterating through
    skymap = butler.get('skyMap',collections='skymaps',dataId={'instrument':'DECam','tract':0,'skymap':'{CLN}_skymap'.format(CLN=cln)})
    
    # this is just a lazy trick for getting the tract-obj from the skymap
    for tract in skymap:
        print(skymap.config)
    
    # I'm writing this such that we can, in principle, change the num of patches
    # since it's not explicitly fixed at being 12x12 (altho it almost always should be)
    xIndex = tract.getNumPatches()[0]
    yIndex = tract.getNumPatches()[1]
    numPatches = xIndex*yIndex
    
    for shear_type in ["noshear", "1p", "2p", "1m", "2m"]:
        
        table_array = []
        
        # load data from each patch (if it exists)
        for patch in range(numPatches):
            
            compute_magnitudes = ['u_inst_psf_flux','u_inst_cmodel_flux','g_inst_psf_flux','g_inst_cmodel_flux','r_inst_psf_flux','r_inst_cmodel_flux','i_inst_psf_flux','i_inst_cmodel_flux','z_inst_psf_flux','z_inst_cmodel_flux']
            tab = export_patch_data(butler,
                                    patch,
                                    flags,
                                    columns,
                                    cln=cln,
                                    compute_magnitudes=compute_magnitudes,
                                    compute_shapes=['sdss','sdss_psf','i_sdss_psf','hsm','hsm_psf','simple'],
                                    dataset_type=f'deep_{shear_type}_Coadd_obj',
                                    collections=f'DECam/processing/meta_4b_{shear_type}',
                                   )
            
            table_array.append(tab)
        
        # and finally write the data to disk    
        data_out = vstack(table_array)
        
        # before exporting, propagate measurements in sdss_err to the sdss_shear
        shape_type = 'sdss'
        xx = data_out[shape_type + '_xx']
        xx_err = data_out[shape_type + '_xx_err']
        yy = data_out[shape_type + '_yy']
        yy_err = data_out[shape_type + '_yy_err']
        xy = data_out[shape_type + '_xy']
        xy_err = data_out[shape_type + '_xy_err']
        
        T = xx + yy
        T_err = np.sqrt(xx_err**2 + yy_err**2)
        
        # for division relative errors add in quadrature, T_err = err_(xx - yy) and err_(2xy)/(2xy) = err_(xy)/xy
        data_out[shape_type + '_e1_err'] = data_out[shape_type + '_e1'] * np.sqrt( (T_err/T)**2 + (T_err/(xx - yy))**2 )
        data_out[shape_type + '_e2_err'] = data_out[shape_type + '_e2'] * np.sqrt( (T_err/T)**2 + (xy_err/xy)**2 )
        
        # also compute the resolvedness here
        data_out[shape_type + '_res'] = 1 - ( ( data_out[shape_type + '_psf_xx'] + data_out[shape_type + '_psf_yy'] )/( data_out[shape_type + '_xx'] + data_out[shape_type + '_yy'] ) )
                
        data_out.write(f'metadetect_export/{cln}_metadetect_{shear_type}.csv', format="ascii.csv", overwrite=True)
        
        
        



