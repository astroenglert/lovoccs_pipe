import sys

import numpy as np
import matplotlib.pyplot as pl

from astropy.wcs import WCS
from astropy import units as u

from lsst.daf.butler import Butler

# some arguments/main statement with functions below as alwys
# this script draws a handy little guide where each patch is labelled on the fov of the cluster
# with the lines of ra/dec separating each patch drawn on

def get_patch_info(tract):
    '''
    Helper function to load information about a tract
    
    Args:
      tract: int; integer specifying the tract-id
    
    Returns:
      id_array: array; array of patch id-numbers
      bbox_dict: dict; dictionary of bbox info keyed by patch-index
    
    '''
    # ususally the skymap should be 12x12, but I'm leaving this general in case something changes
    xMax = tract.getNumPatches()[0]
    yMax = tract.getNumPatches()[1]
    
    # adding the "patch-id" to an array which describes the patches location
    id_array = np.ones((xMax,yMax))
    for i in range(xMax):
        for j in range(yMax):
            id_array[i,j] = tract.getSequentialPatchIndexFromPair((i,j))
    
    # create an array containing the bbox objects for each patch
    bbox_dict = {}
    
    for i in range(xMax):
        for j in range(yMax):
            info = tract.getPatchInfo((i,j))
            bbox_dict[(i,j)] = info.getInnerBBox()

    return id_array,bbox_dict

def draw_patches(ax,id_array,bbox_dict,textparams=None):
    '''
    Helper function to render the patches on-sky (roughly...)
    
    Args:
      ax: matplotlib Axes; axes to draw patches over
      id_array: array; array of patch-ID's
      bbox_dict: dict; dictionary of bbox info keyed by patch index
      textparams: dict; dictionary of parameter to pass to ax.annotate, defaults to None
    
    Returns:
      None

    '''
    
    # first just a quick lazy way to get the number of patches
    xMax = len(id_array[:,0])
    yMax = len(id_array[0,:])
    
    # now let's iterate through each patch and draw the patch boundaries on the axis
    corners = []
    for i in range(xMax):
        for j in range(yMax):
            bbox = bbox_dict[(i,j)]
            corners.append([bbox.getMinX(),bbox.getMinY()])
            center = bbox.getCenter()
            ax.annotate(f'({i},{j})',(center[0],center[1]),verticalalignment='top',**textparams)
            ax.annotate(int(id_array[i,j]),(center[0],center[1]),verticalalignment='bottom',color='red',**textparams)
    bbox_edge = bbox_dict[xMax-1,yMax-1]
    corners.append([bbox_edge.getMaxX(),bbox_edge.getMaxY()])
    
    # instead of drawing line, just let the numpy grid take care of it for you!
    corners = np.array(corners)
    x_ticks = np.unique(corners[:,0])
    y_ticks = np.unique(corners[:,1])
    ticks_deg = wcs.all_pix2world(x_ticks,y_ticks,0)
    ax.coords[0].set_ticks(ticks_deg[0]*u.deg)
    ax.coords[1].set_ticks(ticks_deg[1]*u.deg)
    ax.coords.grid(color='gray',alpha=0.8,linestyle='--')
    
    # just generic formatting here, getting units of degrees with one decimal
    ax.coords[0].set_major_formatter('d.d')
    ax.coords[1].set_major_formatter('d.d')
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('DEC')
    
    ax.set_xlim((min(x_ticks),max(x_ticks)))
    ax.set_ylim((min(y_ticks),max(y_ticks)))
    
# run if not initialized from import (I'll be using these functions a lot, so this will be useful)
if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise Exception("Proper Usage: python draw_skymap.py REPO CLN")
    
    repo_name = sys.argv[1]
    cln = sys.argv[2]
    
    butler = Butler(repo_name)

    # retrieve the skymap from the butler
    skymap = butler.get('skyMap',collections='skymaps',dataId={'instrument':'DECam','tract':0,'skymap':'{CLN}_skymap'.format(CLN=cln)})
    
    # a lazy trick to get the 'tract' object, which contains info about patches
    for tract in skymap:
        print(skymap.config)

    # creating the wcs
    wcs = WCS(tract.getWcs().getFitsMetadata())
    ax = pl.subplot(projection=wcs)
    
    # getting patch info and drawing
    id_array,bbox_array = get_patch_info(tract)
    draw_patches(ax,id_array,bbox_array,textparams={'horizontalalignment':'center','fontsize':'x-small'})
    pl.savefig("skymap_indices.png",dpi=480,bbox_inches='tight')
    
    #TODO query aladin/simbad/somewhere to get an sdss or decals image to put in the background
    # this actually won't be too helpful since the fov is much larger than DECam's footprint


    
    




