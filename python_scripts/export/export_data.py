import os
import time
import math
import glob
import sys
import json

from importlib import resources as impresources
from pathlib import Path

import numpy as np

import astropy as ap
import astropy.io.fits as pyfits
import astropy.units as u

from astropy.wcs import WCS

from astropy.visualization import ManualInterval
from astropy.visualization import AsinhStretch
from astropy.visualization import LogStretch
from astropy.visualization import make_lupton_rgb

from astropy.table import Table, vstack

from astropy.coordinates import SkyCoord

import pandas as pd
import matplotlib.pyplot as pl

from lsst.daf.butler import Butler

# homebrew modules below
from . import export_config
export_config = impresources.files(export_config)


def export_patch_data(butler,patch,flags,columns,cln='A85',compute_magnitudes=['r_inst_psf_flux','r_inst_cmodel_flux'],compute_shapes=['sdss','sdss_psf','hsm','hsm_psf','i_sdss_psf'],dataset_type='deepCoadd_obj',collections='DECam/processing/coadd_3c'):
    '''
    A function for exporting data from a specified patch, given a dictionary of flags and columns. The formatting of the flag dictionaries is {"BAND:TABLE_NAME:COLUMN_NAME":"TRUE/FALSE"} and the column dictionaries are structured similarly" {"BAND:TABLE_NAME:COLUMN_NAME":"OUTPUT_COLUMN_NAME"}. This was done to make exporting custom measurements a bit easier (you're welcome Soren!).
    
    Args:
        butler: LSST Butler; a butler for the repository you're exporting from
        patch: int; an index specifying the patch-index to export from
        flags: dict; a dictionary of flags used for filtering objects
        columns: dict; a dictionary of columns to export
        flux_types: array; a list of output fluxes to transform into magnitudes
        compute_shapes: array; a list of moment-types to compute shapes for (xx,yy,xy -> e1/e2)
    
    Returns:
        patch_table: Astropy Table; a table containing the specified columns
    
    '''
    
    patch_table = Table()
    print(f'Trying to export Patch {patch}')
    # if the deepCoaddObj table doesn't exist, or the flags are missing, skip the patch!
    try:
        deepCoaddObj = butler.get(dataset_type,dataId={'instrument':'DECam','skymap':'{CLN}_skymap'.format(CLN=cln),'tract':0,'patch':patch},collections=collections)

        # collect the object flags
        select_objects = np.ones(len(deepCoaddObj)).astype(bool)
        for input_col,output_col in flags.items():
            input_col = input_col.split(":")
            if output_col == "False":
                select_objects &= ~deepCoaddObj[input_col[1]][input_col[0]][input_col[2]]
            else:
                select_objects &= deepCoaddObj[input_col[1]][input_col[0]][input_col[2]]
    
    except:
        print(f'Failed to export Patch {patch}, skipping...')
        # if the above fails, export an empty table
        return patch_table
    
    # now select objects from the table and begin exporting
    deepCoaddObj_selected = deepCoaddObj[select_objects]
    table_size = len(deepCoaddObj_selected)
    
    # try exporting each column from deepCoaddObj, if it fails export an array of NaN
    for input_col,output_col in columns.items():
    
        input_col = input_col.split(":")
        try:
            values = deepCoaddObj_selected[input_col[1]][input_col[0]][input_col[2]]
        except:
            print('Unable to export %s from %s in band %s... assigning NaN'%(input_col[2],input_col[1],input_col[0]))
            values = np.full(table_size,np.nan)
        
        patch_table[output_col] = values
        
        #TODO should we separate these somewhere in the json? Or is this a sufficient way of identifying magnitudes?
        if output_col in compute_magnitudes:
            
            # try to load photocalib for the band
            try:
                CoaddPhotoCalib = butler.get('deepCoadd_calexp.photoCalib',dataId={'instrument':'DECam','skymap':'{CLN}_skymap'.format(CLN=cln),'tract':0,'patch':patch,'band':input_col[0]},collections=collections)
                mag = CoaddPhotoCalib.instFluxToMagnitudeArray(values.values,x=patch_table['x'].data,y=patch_table['y'].data) 
                
                # load the fluxerr for this band
                #TODO is the standard to always add 'Err' to make error columns? If not this will need to be changed
                fluxerr = deepCoaddObj_selected[input_col[1]][input_col[0]][input_col[2] + 'Err']
                magerr = ( 2.5/np.log(10) ) * (fluxerr/values)
            
            except:
                print('Unable to load the photocalib in band %s for patch %s'%(input_col[0],patch))
                fluxerr = np.full(table_size,np.nan)
                mag = np.full(table_size,np.nan)
                magerr = np.full(table_size,np.nan)
            
            #TODO potential minor bug, this is only assigned if the magnitude is computed, should always be assigned
            patch_table[output_col + 'err'] = fluxerr
            
            #TODO this is a lazy way of getting the flux-type, will need to change in the future
            flux_type = output_col.split('_')[2]
            patch_table[input_col[0] + '_' + flux_type + '_mag'] = mag
            patch_table[input_col[0] + '_' + flux_type + '_magerr'] = magerr
    
    # finally, compute shapes
    for shape in compute_shapes:
        
        xx = patch_table[shape + '_xx']
        yy = patch_table[shape + '_yy']
        xy = patch_table[shape + '_xy']
        
        patch_table[shape + '_e1'] = (xx - yy)/(xx + yy)
        patch_table[shape + '_e2'] = (2*xy)/(xx + yy)
    
    # change the units of the ra/dec columns to degress
    patch_table['ra'] = patch_table['ra'] * 180/np.pi
    patch_table['dec'] = patch_table['dec'] * 180/np.pi
    
    # add a unique identifier
    #TODO should we do LVS or LV[XRAY RANK NUM]; e.g. for A85 LVS 0039... v. LV3 0039...
    coords = SkyCoord(ra=patch_table['ra'],dec=patch_table['dec'],unit='deg')
    ra_str = coords.ra.to_string(u.hourangle,precision=2,sep='',pad=True)
    dec_str = coords.dec.to_string(u.degree,precision=2,sep='',pad=True,alwayssign=True)
    
    object_id = np.char.add('LVS J',ra_str)
    object_id = np.char.add(object_id,dec_str)
    patch_table.add_column(object_id,name='ID',index=0)
    
    return patch_table


# by default we clip and stretch the data before displaying it
def display(mat, vmin, vmax, tag):
    interval = ManualInterval(vmin, vmax)
    mat_1 = interval(mat)
    stretch = AsinhStretch(0.0004)
    mat_2 = stretch(mat_1)

    return mat_2


# define a function here to load and save fits-files of a given patch range
# will export the deepCoadd and deepCoadd_calexp (deepcoadd+final round bkg-subtraction)
def export_patches(band,x_list,y_list,patch_width=4000,gap_width=100,dataset_type='deepCoadd',wcs_patch=[5,5],collection='DECam/processing/coadd_3a',tag=None):
    '''
    Helper function to export patches of dataset_type from the repo
    
    Args:
      band: string; the band for the dataset
      x_list: array; an array specifying the x-indices of the patch
      y_list: array; an array specifying the y-indices of the patch
      patch_width: int; dimensions of the patch in DECam pixels, defaults to '4000'
      gap_width: int; the width of the outer-bbox around each patch, defaults to '100'
      dataset_type: string; dataset_type to export, defaults to 'deepCoadd'
      wcs_patch: array; array specifying the patch to use as reference for the WCS, defaults to [5,5]
      collection: string; collection in repo to load from, defaults to 'DECam/processing/coadd_3a
      tag: string; tag to append to the end of the filename, defaults to None
    
    Returns:
      None
    
    '''
    patchIndices = {}
    image_whole = np.zeros((len(x_list)*patch_width, len(y_list)*patch_width))
    # first let's get a dictionary of the patches in a square of x_list x y_list
    for i in x_list:
        for j in y_list:
            patchIndices[tract.getSequentialPatchIndexFromPair((i,j))]=(i,j)

    # try to load the patch and continue if it doesn't exist
    for patch in list(patchIndices):
        try:
            image_path = butler.getURI(dataset_type,collections=collection,dataId={'instrument':'DECam','band':band,'skymap':'{CLN}_skymap'.format(CLN=cln),'patch':patch,'tract':0}).ospath
        except:
            print("Patch {PTC} doesn't exist for {BAND}, skipping...".format(PTC=patchIndices[patch],BAND=band))
            print("Zero will be assigned to the img here!")
            continue
            
        row_ind = patchIndices[patch][0] - min(x_list)
        col_ind = patchIndices[patch][1] - min(y_list)
        
        x_ind = patchIndices[patch][0]
        y_ind = patchIndices[patch][1]
        
        image = pyfits.getdata(image_path, 1)
        if x_ind==0 and y_ind==0:
            image = image[gap_width*0:patch_width+gap_width*0, gap_width*0:patch_width+gap_width*0]
        elif x_ind==0 and y_ind!=0:
            image = image[gap_width:patch_width+gap_width, gap_width*0:patch_width+gap_width*0]
        elif x_ind!=0 and y_ind==0:
             image = image[gap_width*0:patch_width+gap_width*0, gap_width:patch_width+gap_width]
        else:
            image = image[gap_width:patch_width+gap_width, gap_width:patch_width+gap_width]
     
        print('Pasting to whole image...')
        image_whole[(row_ind*patch_width):((row_ind+1)*patch_width), (col_ind*patch_width):((col_ind+1)*patch_width)] = np.transpose(image)
         
    # now that we have everything saved in an array, we just need to create the correct wcs
    # by default this is done from patch 55
    pmin_x=min(x_list)
    pmin_y=min(y_list)

    x_wcs = wcs_patch[0]
    y_wcs = wcs_patch[1]
    
    patch = tract.getSequentialPatchIndexFromPair((x_wcs,y_wcs)) 
    try:
        image_path = butler.getURI(dataset_type,collections=collection,dataId={'instrument':'DECam','band':band,'skymap':'{CLN}_skymap'.format(CLN=cln),'patch':patch,'tract':0}).ospath
    except:
        print("Patch {PTC} doesn't exist for {BAND}, so I'll skip making the fits for this band...".format(PTC=patchIndices[patch],BAND=band))
        return False

    hdul = pyfits.open(image_path)
    hdu = hdul[1]
    w_old = WCS(hdu.header)
    hdul.close()

    # updating the reference pixel and copying any other data we need
    # let's always use the WCS from the center-patch, which always exists!
    w_new = WCS(naxis=2)
    w_new.wcs.crpix = [
                        w_old.wcs.crpix[0]-100+(x_wcs-pmin_x)*4000, 
                        w_old.wcs.crpix[1]-100+(y_wcs-pmin_y)*4000,
                    ]
    w_new.wcs.crval = w_old.wcs.crval
    w_new.wcs.ctype = w_old.wcs.ctype
    w_new.wcs.cd = w_old.wcs.cd
    w_new.wcs.mjdobs = w_old.wcs.mjdobs
    w_new.wcs.dateobs = w_old.wcs.dateobs
    w_new.wcs.radesys = w_old.wcs.radesys

    # And now we can finally write the image
    print('===writing %s band whole image to output...==='%band)    
    header = w_new.to_header()
    hdu = pyfits.PrimaryHDU(image_whole.T, header=header)
    if tag==None:
        output_str = "combine_patch_color_output/%s_%s%d%d-%d%d_%s.fits"%(
                        cln, 
                        band, 
                        min(x_list), 
                        min(y_list), 
                        max(x_list), 
                        max(y_list),
                        dataset_type,
                    )
    else:
        output_str = "combine_patch_color_output/%s_%s%d%d-%d%d_%s_%s.fits"%(
                        cln, 
                        band, 
                        min(x_list), 
                        min(y_list), 
                        max(x_list), 
                        max(y_list),
                        dataset_type,
                        tag
                    )
    
    hdu.writeto(output_str,overwrite=True)
    return True


# Let's also use two stretches, one for color-correct (lupton) and another for not color-correct (but nice looking)
def draw_rgb(R,G,B,patches,tag='deepCoadd',lupton=False):
    '''
    Helper function to draw an RGB image of the data
    
    Args:
      R/G/B: string; string specifying the band to load into each channel
      patches: string; the range of patches (e.g. '33-88')
      tag: string; tag to append to the end of the output file
      lupton: bool; enable Lupton coloring?
    
    Returns
      None
    
    '''
    fits_image_filename_1 = "combine_patch_color_output/{CLN}_{R}{PATCHES}_{TAG}.fits".format(CLN=cln,PATCHES=patches,R=R,TAG=tag)
    fits_image_filename_2 = "combine_patch_color_output/{CLN}_{G}{PATCHES}_{TAG}.fits".format(CLN=cln,PATCHES=patches,G=G,TAG=tag)
    fits_image_filename_3 = "combine_patch_color_output/{CLN}_{B}{PATCHES}_{TAG}.fits".format(CLN=cln,PATCHES=patches,B=B,TAG=tag)
    out = "combine_patch_color_output/{CLN}_{PATCHES}_{R}{G}{B}_{TAG}".format(CLN=cln,PATCHES=patches,R=R,G=G,B=B,TAG=tag)
    
    print("Loading FITS image...")
    if os.path.exists(fits_image_filename_1) and os.path.exists(fits_image_filename_2) and os.path.exists(fits_image_filename_3):
        print("Files all exist!")
    else:
        print("Some files could be missing...\nExiting...")
        return False
    data_1 = (pyfits.getdata(fits_image_filename_1, 0))
    data_2 = (pyfits.getdata(fits_image_filename_2, 0))
    data_3 = (pyfits.getdata(fits_image_filename_3, 0))
    
    # for pretty-pictures, -0.1 and 250 seem to work best
    # but for science -0.06 and 190 best showcase our issues with bckg-subtr and colors
    
    # if Lupton coloring is enabled, we use a different transformation
    print("Scaling image...")
    
    if lupton:
        data_all = make_lupton_rgb(data_1,data_2,data_3,minimum=-0.1,Q=8,stretch=1)
        pl.imsave("%s_lupton.png"%out, 
            np.flipud(data_all),
           )
        print("Writing %s.png..."%out)
    else:
        data_1_scaled = display(data_1, -0.1, 250, 'fig1')*255
        data_2_scaled = display(data_2, -0.1, 250, 'fig2')*255
        data_3_scaled = display(data_3, -0.1, 250, 'fig3')*255
        
        data_1_scaled = np.flipud(data_1_scaled)
        data_2_scaled = np.flipud(data_2_scaled)
        data_3_scaled = np.flipud(data_3_scaled)
        data_all = np.dstack((
                                data_1_scaled, 
                                data_2_scaled, 
                                data_3_scaled,
                            )).astype(np.uint8)
        print("Writing %s.png..."%out)
        pl.imsave("%s.png"%out, 
            data_all,
           )
    return True


if __name__ == '__main__':
    
    # checking for correct usage
    if len(sys.argv)!=2:
        print("python this_script.py cln")
        sys.exit(1)
    
    cln = sys.argv[1]
    
    # create a butler to read/write data with
    butler = Butler('repo/repo')
    
    # load dictionaries of flags and columns
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
    
    table_array = []
    for patch in range(numPatches):
        
        compute_magnitudes = ['u_inst_psf_flux','u_inst_cmodel_flux','g_inst_psf_flux','g_inst_cmodel_flux','r_inst_psf_flux','r_inst_cmodel_flux','i_inst_psf_flux','i_inst_cmodel_flux','z_inst_psf_flux','z_inst_cmodel_flux']
        tab = export_patch_data(butler,patch,flags,columns,cln=cln,compute_magnitudes=compute_magnitudes,compute_shapes=['sdss','sdss_psf','hsm','hsm_psf','i_sdss_psf'])
        table_array.append(tab)
        
    data_out = vstack(table_array)
    
    data_out.write("read_catalog_all_output/{CLN}_00-1111_all.csv".format(CLN=cln), format="ascii.csv", overwrite=True)
    
    # by default, let's export fits of 2x2,4x4,6x6, and the 12x12 fov
    # avoid assuming a particular number of patches
    patch_22 = [np.arange(int(xIndex/2)-1,int(xIndex/2)+1),np.arange(int(yIndex/2)-1,int(yIndex/2)+1)]
    ranges_22 = str(min(patch_22[0])) + str(min(patch_22[0])) + '-' + str(max(patch_22[1])) + str(max(patch_22[1]))
    patch_44 = [np.arange(int(xIndex/2)-2,int(xIndex/2)+2),np.arange(int(yIndex/2)-2,int(yIndex/2)+2)]
    ranges_44 = str(min(patch_44[0])) + str(min(patch_44[0])) + '-' + str(max(patch_44[1])) + str(max(patch_44[1]))
    patch_66 = [np.arange(int(xIndex/2)-3,int(xIndex/2)+3),np.arange(int(yIndex/2)-3,int(yIndex/2)+3)]
    ranges_66 = str(min(patch_66[0])) + str(min(patch_66[0])) + '-' + str(max(patch_66[1])) + str(max(patch_66[1]))
    patch_all = [np.arange(0,xIndex),np.arange(0,yIndex)]
    ranges_all = str(min(patch_all[0])) + str(min(patch_all[0])) + '-' + str(max(patch_all[1])) + str(max(patch_all[1]))
    
    for patches in [ patch_22, patch_44, patch_66, patch_all ]:
        for band in ['g', 'r', 'i', 'z','u']:
            
            # export patches for each band and deepCoadd
            export_patches(band,patches[0],patches[1],dataset_type='deepCoadd',collection='DECam/processing/coadd_3a')

    for patches in [ patch_22, patch_44, patch_66, patch_all ]:
        for band in ['g', 'r', 'i', 'z']:
        
            # export patches for each band and deepCoadd
            export_patches(band,patches[0],patches[1],dataset_type='deepCoadd',collection='DECam/processing/skycorr',tag='skycorr')
    
    
    # now to make some color images from these files...
    
    # let F be the function which cuts and stretches the data to be displayed as an RGB-image
    # _lupton images are "color correct" in the sense that, after applying the stretch, the colors are preserved
    # e.g. F(R)/F(G) = R/G
    # but bright cores appear bloated and, artistically, the image lacks color balance
    # these are the "most correct" for science
    
    # default images are "relative color correct" in the sense that, redder objects before stretch remain redder post-stretch
    # e.g. if R1/G1 > R2/G2, then F(R1)/F(G1) > F(R2)/F(G2)
    # suffer from an alternate problem where bright pixels tend to be "whiten" despite having a distinct color
    # these are best for outreach/presentation
    
    draw_rgb('i','r','g',ranges_22,lupton=True)
    draw_rgb('i','r','g',ranges_44,lupton=True)
    draw_rgb('i','r','g',ranges_66,lupton=True)
    draw_rgb('i','r','g',ranges_22,lupton=False)
    draw_rgb('i','r','g',ranges_44,lupton=False)
    draw_rgb('i','r','g',ranges_66,lupton=False)
    
    draw_rgb('i','r','g',ranges_22,lupton=True,tag='deepCoadd_skycorr')
    draw_rgb('i','r','g',ranges_44,lupton=True,tag='deepCoadd_skycorr')
    draw_rgb('i','r','g',ranges_66,lupton=True,tag='deepCoadd_skycorr')
    draw_rgb('i','r','g',ranges_22,lupton=False,tag='deepCoadd_skycorr')
    draw_rgb('i','r','g',ranges_44,lupton=False,tag='deepCoadd_skycorr')
    draw_rgb('i','r','g',ranges_66,lupton=False,tag='deepCoadd_skycorr')
    
    
    # zir is sometimes useful... so create these just-in-case
    #draw_rgb('z','i','r',ranges_22,lupton=False)
    #draw_rgb('z','i','r',ranges_44,lupton=False)
    #draw_rgb('z','i','r',ranges_66,lupton=False)
    
    # and just for fun print the entire fov? Usually runs out of memory and doesn't contain anything useful
    #draw_rgb('i','r','g',ranges_all,lupton=False)
    #draw_rgb('i','r','g',ranges_all,lupton=True)

