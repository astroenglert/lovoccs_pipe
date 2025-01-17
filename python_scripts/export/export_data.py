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

# define a function here to load and save fits-files of a given patch range
# will export the deepCoadd and deepCoadd_calexp (deepcoadd+final round bkg-subtraction)
def export_patches(band,x_list,y_list,patch_width=4000,gap_width=100,dataset_type='deepCoadd',wcs_patch=[5,5],collection='DECam/processing/coadd_3a',tag=None):
   
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


#TODO this should be refractored to export in a smarter way and pull column_list from a config file somewhere
if __name__ == '__main__':
    
    # checking for correct usage
    if len(sys.argv)!=2:
        print("python this_script.py cln")
        sys.exit(1)
    
    cln = sys.argv[1]
    
    # create a butler to read/write data with
    butler = Butler('repo/repo')
    
    # Y band will (unless explicitly enabled) be disabled... but has to be output as NaN here since later steps
    # depend explicitly (for some reason) on having an exact number of columns -_- rather than looking for headers
    band_list = ['u', 'g', 'r', 'i', 'z', 'Y']
    
    # creating a numpy array for storing all the data
    column_list = ["ra", "dec", "x", "y", "e1", "e2", "res", "sigmae", "rkron", "extendedness", "blendedness", "psf_used", "e1_sdss", "e2_sdss", "e1_psf_sdss", "e2_psf_sdss", "e1_hsm", "e2_hsm", "e1_psf_hsm", "e2_psf_hsm", "i_e1_psf_sdss", "i_e2_psf_sdss"]
    
    # by default, gen3 only outputs fluxes
    for band in band_list:
        column_list.append(band + "_psf_mag")
        column_list.append(band + "_psf_magerr")
        column_list.append(band + "_cmodel_mag")
        column_list.append(band + "_cmodel_magerr")

    data = {}
    for column in column_list:
    # Use dictionary: string -> array to store columns of the catalog
        data[column] = np.array([])
    
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
    
    for patch in range(numPatches):
        
        # first I need to load the deepCoadd_obj for the patch
        try:
            # if this doesn't load or r-band doesn't exist, we skip it!
            deepCoaddObj = butler.get('deepCoadd_obj',dataId={'instrument':'DECam','skymap':'{CLN}_skymap'.format(CLN=cln),'tract':0,'patch':patch},collections='DECam/processing/coadd_3c')
            refTable = deepCoaddObj['ref']['r']
            Sources = deepCoaddObj['forced_src']['r']
            Coadd = deepCoaddObj['meas']['r']
    
        except:
            print("No r-band in {PATCH}, skipping...".format(PATCH=patch))
            continue

        isPrimary = refTable['detect_isPrimary']
        Extendedness = Sources["base_ClassificationExtendedness_value"]
        isGoodFlux = ~Sources["modelfit_CModel_flag"]    
        isGoodFlux &= ~Sources["base_PsfFlux_flag"]    
        selected = isPrimary & isGoodFlux
        
        RA = Sources["coord_ra"]/np.pi*180
        DEC = Sources["coord_dec"]/np.pi*180
    
        X = Coadd["base_SdssCentroid_x"]
        Y = Coadd["base_SdssCentroid_y"]
        isGoodPosition = ~Coadd["base_SdssCentroid_flag"]
        selected &= isGoodPosition
        
        E1 = Coadd["ext_shapeHSM_HsmShapeRegauss_e1"]
        E2 = Coadd["ext_shapeHSM_HsmShapeRegauss_e2"]
    
        RES = Coadd["ext_shapeHSM_HsmShapeRegauss_resolution"]
        SIGMAE = Coadd["ext_shapeHSM_HsmShapeRegauss_sigma"]
    
        RKRON = Coadd["ext_photometryKron_KronFlux_radius"]
        
        # 0 for isolated objects, strongly blended -> 1
        Blendedness = Coadd["base_Blendedness_abs"]
        
        #--------------------------------------------
        PSF_USED = Coadd["calib_psf_used"] * 1
        
        SDSS_SHAPE_I11 = Coadd['base_SdssShape_xx']
        SDSS_SHAPE_I12 = Coadd['base_SdssShape_xy']
        SDSS_SHAPE_I22 = Coadd['base_SdssShape_yy']
        E1_SDSS = (SDSS_SHAPE_I11 - SDSS_SHAPE_I22) / (SDSS_SHAPE_I11 + SDSS_SHAPE_I22)
        E2_SDSS = 2. * SDSS_SHAPE_I12 / (SDSS_SHAPE_I11 + SDSS_SHAPE_I22)
        
        #--------------------------------------------
        SDSS_PSF_I11 = Coadd['base_SdssShape_psf_xx']
        SDSS_PSF_I12 = Coadd['base_SdssShape_psf_xy']
        SDSS_PSF_I22 = Coadd['base_SdssShape_psf_yy']
        E1_PSF_SDSS = (SDSS_PSF_I11 - SDSS_PSF_I22) / (SDSS_PSF_I11 + SDSS_PSF_I22)
        E2_PSF_SDSS = 2. * SDSS_PSF_I12 / (SDSS_PSF_I11 + SDSS_PSF_I22)
        
        HSM_SHAPE_I11 = Coadd['ext_shapeHSM_HsmSourceMoments_xx']
        HSM_SHAPE_I12 = Coadd['ext_shapeHSM_HsmSourceMoments_xy']
        HSM_SHAPE_I22 = Coadd['ext_shapeHSM_HsmSourceMoments_yy']
        E1_HSM = (HSM_SHAPE_I11 - HSM_SHAPE_I22) / (HSM_SHAPE_I11 + HSM_SHAPE_I22)
        E2_HSM = 2. * HSM_SHAPE_I12 / (HSM_SHAPE_I11 + HSM_SHAPE_I22)
        
        HSM_PSF_I11 = Coadd['ext_shapeHSM_HsmPsfMoments_xx']
        HSM_PSF_I12 = Coadd['ext_shapeHSM_HsmPsfMoments_xy']
        HSM_PSF_I22 = Coadd['ext_shapeHSM_HsmPsfMoments_yy']
        E1_PSF_HSM = (HSM_PSF_I11 - HSM_PSF_I22) / (HSM_PSF_I11 + HSM_PSF_I22)
        E2_PSF_HSM = 2. * HSM_PSF_I12 / (HSM_PSF_I11 + HSM_PSF_I22)
        
        try:
            
            # try to collect the i-band shapes for shear calibration
            iCoadd = deepCoaddObj['meas']['i']
            I_SDSS_PSF_I11 = iCoadd['base_SdssShape_psf_xx']
            I_SDSS_PSF_I12 = iCoadd['base_SdssShape_psf_xy']
            I_SDSS_PSF_I22 = iCoadd['base_SdssShape_psf_yy']
            I_E1_PSF_SDSS = (I_SDSS_PSF_I11 - I_SDSS_PSF_I22) / (I_SDSS_PSF_I11 + I_SDSS_PSF_I22)
            I_E2_PSF_SDSS = 2. * I_SDSS_PSF_I12 / (I_SDSS_PSF_I11 + I_SDSS_PSF_I22)
        
        except:
            
            print("No i-band in {PATCH}, setting psf-shape to NaN...".format(PATCH=patch))
            I_E1_PSF_SDSS = np.full_like(selected, np.nan, dtype=np.double)
            I_E2_PSF_SDSS = np.full_like(selected, np.nan, dtype=np.double)

        # now collecting the photometry
        for band in band_list:
            
            try:
                Sources = deepCoaddObj['forced_src'][band]
                CoaddPhotoCalib = butler.get('deepCoadd_calexp.photoCalib',dataId={'instrument':'DECam','skymap':'{CLN}_skymap'.format(CLN=cln),'tract':0,'patch':patch,'band':band},collections='DECam/processing/coadd_3c')
            except:
                print("No photocalib or source catalog for {BAND} in {PATCH}, setting to NaN...".format(BAND=band,PATCH=patch))
                data[band + "_psf_mag"] = np.append( data[band + "_psf_mag"], np.full_like(selected[selected], np.nan, dtype=np.double) )
                data[band + "_psf_magerr"] = np.append( data[band + "_psf_magerr"], np.full_like(selected[selected], np.nan, dtype=np.double) )
                data[band + "_cmodel_mag"] = np.append( data[band + "_cmodel_mag"], np.full_like(selected[selected], np.nan, dtype=np.double) )
                data[band + "_cmodel_magerr"] = np.append( data[band + "_cmodel_magerr"], np.full_like(selected[selected], np.nan, dtype=np.double) )
                continue
            
            #TODO this is all correct... but omits error in the scale
            # shouldn't be an issue (scale is fixed by mzp=27 and measurement error is dom noise)
            # but eventually we should fix this
            Psf_Mags = CoaddPhotoCalib.instFluxToMagnitudeArray(Sources["slot_PsfFlux_instFlux"].values,x=Sources["base_SdssCentroid_x"].values,y=Sources["base_SdssCentroid_y"].values) 
            CModel_Mags = CoaddPhotoCalib.instFluxToMagnitudeArray(Sources["modelfit_CModel_instFlux"].values,x=Sources["base_SdssCentroid_x"].values,y=Sources["base_SdssCentroid_y"].values) 
            
            data[band + "_psf_mag"] = np.append( data[band + "_psf_mag"], Psf_Mags[selected].value )
            data[band + "_psf_magerr"] = np.append( data[band + "_psf_magerr"], ( 2.5/np.log(10) ) * Sources[selected]["slot_PsfFlux_instFluxErr"]/Sources[selected]["slot_PsfFlux_instFlux"] )
            data[band + "_cmodel_mag"] = np.append( data[band + "_cmodel_mag"], CModel_Mags[selected].value )
            data[band + "_cmodel_magerr"] = np.append( data[band + "_cmodel_magerr"], ( 2.5/np.log(10) ) * Sources[selected]["modelfit_CModel_instFluxErr"]/Sources[selected]["modelfit_CModel_instFlux"] )
            
        data["ra"] = np.append(data["ra"], RA[selected])
        data["dec"] = np.append(data["dec"], DEC[selected])
        
        data['x'] = np.append(data['x'], X[selected])
        data['y'] = np.append(data['y'], Y[selected])
        data["e1"] = np.append(data["e1"], E1[selected])
        data["e2"] = np.append(data["e2"], E2[selected])
        data["res"] = np.append(data["res"], RES[selected])
        data["sigmae"] = np.append(data["sigmae"], SIGMAE[selected])
        data["rkron"] = np.append(data["rkron"], RKRON[selected])
        
        data["extendedness"] = np.append(data["extendedness"], Extendedness[selected])
        data["blendedness"] = np.append(data["blendedness"], Blendedness[selected])
        
        data["psf_used"] = np.append(data["psf_used"], PSF_USED[selected])
        
        data["e1_sdss"] = np.append(data["e1_sdss"], E1_SDSS[selected])
        data["e2_sdss"] = np.append(data["e2_sdss"], E2_SDSS[selected])
        
        data["e1_psf_sdss"] = np.append(data["e1_psf_sdss"], E1_PSF_SDSS[selected])
        data["e2_psf_sdss"] = np.append(data["e2_psf_sdss"], E2_PSF_SDSS[selected])
        
        data["e1_hsm"] = np.append(data["e1_hsm"], E1_HSM[selected])
        data["e2_hsm"] = np.append(data["e2_hsm"], E2_HSM[selected])
        
        data["e1_psf_hsm"] = np.append(data["e1_psf_hsm"], E1_PSF_HSM[selected])
        data["e2_psf_hsm"] = np.append(data["e2_psf_hsm"], E2_PSF_HSM[selected])
        
        data["i_e1_psf_sdss"] = np.append(data["i_e1_psf_sdss"], I_E1_PSF_SDSS[selected])
        data["i_e2_psf_sdss"] = np.append(data["i_e2_psf_sdss"], I_E2_PSF_SDSS[selected])
    
    print("\n\n" + '#'*30 + "\nBuilding output table...")
    
    data_out = Table()
    
    for column in column_list:
        print("Column %s: len: %s"%(column, len(data[column]) ) )
        
        data_out[column] = data[column]
    
    # need to overhaul this with a better system eventually...
    # Self-defined id (starting from 1)
    data_out["idn"] = np.array([i for i in range(1, len(data_out)+1) ])
    
    data_out.write("read_catalog_all_output/{CLN}_00-1111_all.csv".format(CLN=cln), format="ascii.csv", overwrite=True)
    
    # That takes care of the old "read_catalog_all" step from our old pipeline
    # read merge full will not be carred over since those catalogs can be easily queried by the butler

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
    
    for patches in [ patch_22, patch_44, patch_66 ]:
        for band in ['g', 'r', 'i', 'z','u']:
            
            # export patches for each band and deepCoadd
            export_patches(band,patches[0],patches[1],dataset_type='deepCoadd',collection='DECam/processing/coadd_3a')

    for patches in [ patch_22, patch_44, patch_66 ]:
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

