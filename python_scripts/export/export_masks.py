# helper function to create either a fits or png for the given patch range
import os
import time
import math
import glob
import sys

import numpy as np

import astropy as ap
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.visualization import ManualInterval
from astropy.visualization import AsinhStretch
from astropy.visualization import LogStretch
from astropy.visualization import make_lupton_rgb
from astropy.table import Table

import pandas as pd
import matplotlib.pyplot as pl

from lsst.daf.butler import Butler

# by default we clip and stretch the data before displaying it
def display(mat, vmin, vmax, tag):
    interval = ManualInterval(vmin, vmax)
    mat_1 = interval(mat)
    stretch = AsinhStretch(0.0004)
    mat_2 = stretch(mat_1)

    return mat_2

def collect_masks(butler,tract,cln,band,x_list,y_list,dataset_type='deepCoadd',patch_width=4000,gap_width=100,mask_planes=['BAD','DETECTED','DETECTED_NEGATIVE','SAT'],wcs_patch=[5,5],collection='DECam/processing/coadd_3a',tag=None):
    '''
    Function to extract a mask from filepath located at hdul_index.
    Will extract all masks specified by the bits above, by default for coadds it will extrct these useful ones:
    0 -> MP_BAD
    1 -> MP_SAT
    5 -> MP_DETECTED
    6 -> MP_DETECTED_NEGATIVE
    
    For more masks, print the header of hdul_index.
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
            image = butler.get(dataset_type,collections=collection,dataId={'instrument':'DECam','band':band,'skymap':'{CLN}_skymap'.format(CLN=cln),'patch':patch,'tract':0})
        except:
            print("Patch {PTC} doesn't exist for {BAND}, skipping...".format(PTC=patchIndices[patch],BAND=band))
            print("Zero will be assigned to the img here!")
            continue
            
        row_ind = patchIndices[patch][0] - min(x_list)
        col_ind = patchIndices[patch][1] - min(y_list)
        
        x_ind = patchIndices[patch][0]
        y_ind = patchIndices[patch][1]
        
        # now load the mask and select the relevant planes
        masks = image.getMask()
        mask_plane_dict = masks.getMaskPlaneDict()
        image = None
        for plane in mask_planes:
            bit = masks.getMaskPlane(plane)
            if image is None:
                image = masks.array >> bit & 1
            else:
                new_mask = masks.array >> bit & 1
                image = new_mask | image
        
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
    
    if not os.path.exists('masks/'):
        os.makedirs('masks/')
    if tag==None:
        output_str = "masks/%s_%s%d%d-%d%d_%s.fits"%(
                        cln, 
                        band, 
                        min(x_list), 
                        min(y_list), 
                        max(x_list), 
                        max(y_list),
                        dataset_type,
                    )
    else:
        output_str = "masks/%s_%s%d%d-%d%d_%s_%s.fits"%(
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
    
    return None


# define a function here to load and save fits-files of a given patch range
# will export the deepCoadd and deepCoadd_calexp (deepcoadd+final round bkg-subtraction)
def export_patches(butler,tract,cln,band,x_list,y_list,patch_width=4000,gap_width=100,dataset_type='deepCoadd',wcs_patch=[5,5],collection='DECam/processing/coadd_3a',tag=None):
   
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
    
    print("Loading FITS images...")
    data_1 = (pyfits.getdata(R, 0))
    data_2 = (pyfits.getdata(G, 0))
    data_3 = (pyfits.getdata(B, 0))
    
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


if __name__=='__main__':
    
    print(len(sys.argv))
    # two usage cases: first export specific fits patches and second is export specific color image
    if not (len(sys.argv)==7 or len(sys.argv)==8):  
        print("Improper usage, try either: ")
        print("python export_patches.py cln repo coll blp trp band OR python export_patches.py cln repo coll blp trp band dataset_type")
        print("If dataset_type is not specified, it defaults to deepCoadd (really most people won't need to use that option...")
        sys.exit(1)
    
    cln = sys.argv[1]
    repo = sys.argv[2]
    coll = sys.argv[3]
    blp = sys.argv[4]
    trp = sys.argv[5]
    band = sys.argv[6]
    
    if len(sys.argv)==7:
        dataset_type='deepCoadd'
    else:
        dataset_type = sys.argv[7]
    
    # create a butler to read/write data with
    butler = Butler(repo)

    # load the skymap
    # for run_steps, we'll always save this in
    skymap = butler.get('skyMap',collections='skymaps',dataId={'instrument':'DECam','tract':0,'skymap':'{CLN}_skymap'.format(CLN=cln)})

    # this is just a lazy trick for getting the tract-obj from the skymap
    for tract in skymap:
        print(skymap.config)

    # I'm writing this such that we can, in principle, change the num of patches
    # since it's not explicitly fixed at being 12x12 (altho it almost always should be)
    xIndex = tract.getNumPatches()[0]
    yIndex = tract.getNumPatches()[1]
    numPatches = xIndex*yIndex

    # input indices separated by string
    xIndexMin = blp.split(',')[0]
    xIndexMax = trp.split(',')[0]
    yIndexMin = blp.split(',')[1]
    yIndexMax = trp.split(',')[1]

    patches = [np.arange(int(xIndexMin),int(xIndexMax)+1),np.arange(int(yIndexMin),int(yIndexMax)+1)]
    print(patches)
    ranges = str(min(patches[0])) + str(min(patches[0])) + '-' + str(max(patches[1])) + str(max(patches[1]))
    
    #export_patches(butler,tract,cln,band,patches[0],patches[1],dataset_type=dataset_type,wcs_patch=[5,5],collection=coll)
    collect_masks(butler,tract,cln,band,patches[0],patches[1],dataset_type=dataset_type,wcs_patch=[5,5],collection=coll)

