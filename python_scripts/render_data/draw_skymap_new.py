import sys

import numpy as np
import matplotlib.pyplot as pl
from astroquery.hips2fits import hips2fits
from matplotlib.colors import Colormap

from astropy.wcs import WCS
from astropy import units as u
from astroquery.ned import Ned
from astropy.visualization.wcsaxes import SphericalCircle

from lsst.daf.butler import Butler

# some arguments/main statement with functions below as alwys
# this script draws a handy little guide where each patch is labelled on the fov of the cluster
# with the lines of ra/dec separating each patch drawn on

# function to query hips for a cutout of CLN
def get_cluster_cutout(cln="A85",catalog="CDS/P/DESI-Legacy-Surveys/DR10/r",dims=2000,fov=4,hipsQuery={}):
    '''
    Helper function to query HiPS for a cutout of the cluster
    
    Args:
      cln: string; cluster name
      catalog: string; catalog to load data from, defaults to CDS/P/DESI-Legacy-Surveys/DR10/r
      dims: int; width of the .fits image
      fov: int; field of view of the .fits image in degrees
      hipsQuery: dict; arguments to pass to the hipsQuery, defaults to {}
    
    
    '''
    # query for object coordinates
    out = Ned.query_object(cln)
    ra = out[0]['RA']
    dec = out[0]['DEC']
    
    # create the wcs
    w = WCS(header={
    'NAXIS1': dims,         # Width of the output fits/image
    'NAXIS2': dims,         # Height of the output fits/image
    'WCSAXES': 2,           # Number of coordinate axes
    'CRPIX1': dims/2,       # Pixel coordinate of reference point
    'CRPIX2': dims/2,       # Pixel coordinate of reference point
    'CDELT1': -fov/dims,    # [deg] Coordinate increment at reference point
    'CDELT2': -fov/dims,    # [deg] Coordinate increment at reference point
    'CUNIT1': 'deg',        # Units of coordinate increment and value
    'CUNIT2': 'deg',        # Units of coordinate increment and value
    'CTYPE1': 'RA---TAN',   # ra, tangent projection (LSP default)
    'CTYPE2': 'DEC--TAN',   # dec, tangent projection (LSP default)
    'CRVAL1': ra,           # [deg] Coordinate value at reference point
    'CRVAL2': dec,          # [deg] Coordinate value at reference point
    })
    
    # query hips for the cutout
    result = hips2fits.query_with_wcs(
    hips=catalog,
    wcs=w,
    get_query_payload=False,
    format='jpg',
    min_cut=0.5,
    max_cut=99.91,
    cmap=Colormap('gray'),
    stretch='log',
    **hipsQuery,
    )
    
    return w,result

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
            txt1 = ax.annotate(f'({i},{j})',(center[0],center[1]),verticalalignment='top',**textparams)
            txt2 = ax.annotate(int(id_array[i,j]),(center[0],center[1]),verticalalignment='bottom',color='red',**textparams)
            txt1.set_alpha(0.5)
            txt2.set_alpha(0.5)
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
    
    
decam_rds = [
[ [ -0.45920, 0.97905],[-0.45942, 0.82986],[-0.16033, 0.82926],[-0.16019, 0.97832],[-0.45920, 0.97905 ], ],
    [ [ -0.14890, 0.97826],[-0.14931, 0.82920],[0.14998, 0.82855],[0.15021, 0.97761],[-0.14890, 0.97826 ], ],
    [ [ 0.16164, 0.97830],[0.16166, 0.82924],[0.46074, 0.82956],[0.46066, 0.97876],[0.16164, 0.97830 ], ],
    [ [ -0.61452, 0.81569],[-0.61483, 0.66644],[-0.31585, 0.66555],[-0.31557, 0.81475],[-0.61452, 0.81569 ], ],
    [ [ -0.30426, 0.81477],[-0.30445, 0.66560],[-0.00501, 0.66538],[-0.00500, 0.81456],[-0.30426, 0.81477 ], ],
    [ [ 0.00637, 0.81453],[0.00626, 0.66535],[0.30570, 0.66534],[0.30561, 0.81452],[0.00637, 0.81453 ], ],
    [ [ 0.31682, 0.81455],[0.31697, 0.66537],[0.61596, 0.66602],[0.61574, 0.81527],[0.31682, 0.81455 ], ],
    [ [ -0.76975, 0.65225],[-0.77000, 0.50289],[-0.47127, 0.50175],[-0.47103, 0.65112],[-0.76975, 0.65225 ], ],
    [ [ -0.45970, 0.65125],[-0.46009, 0.50189],[-0.16069, 0.50110],[-0.16046, 0.65052],[-0.45970, 0.65125 ], ],
    [ [ -0.14919, 0.65049],[-0.14933, 0.50108],[0.15039, 0.50099],[0.15035, 0.65039],[-0.14919, 0.65049 ], ],
    [ [ 0.16160, 0.65050],[0.16161, 0.50108],[0.46103, 0.50141],[0.46086, 0.65077],[0.16160, 0.65050 ], ],
    [ [ 0.47220, 0.65078],[0.47238, 0.50140],[0.77115, 0.50233],[0.77097, 0.65169],[0.47220, 0.65078 ], ],
    [ [ -0.92499, 0.48907],[-0.92522, 0.33964],[-0.62677, 0.33783],[-0.62642, 0.48735],[-0.92499, 0.48907 ], ],
    [ [ -0.61531, 0.48753],[-0.61562, 0.33802],[-0.31646, 0.33697],[-0.31624, 0.48658],[-0.61531, 0.48753 ], ],
    [ [ -0.30483, 0.48667],[-0.30502, 0.33703],[-0.00523, 0.33667],[-0.00518, 0.48633],[-0.30483, 0.48667 ], ],
    [ [ 0.00615, 0.48608],[0.00613, 0.33640],[0.30592, 0.33665],[0.30580, 0.48628],[0.00615, 0.48608 ], ],
    [ [ 0.31719, 0.48648],[0.31722, 0.33685],[0.61639, 0.33753],[0.61628, 0.48705],[0.31719, 0.48648 ], ],
    [ [ 0.62771, 0.48699],[0.62793, 0.33748],[0.92636, 0.33915],[0.92622, 0.48858],[0.62771, 0.48699 ], ],
    [ [ -0.92511, 0.32506],[-0.92527, 0.17559],[-0.62687, 0.17365],[-0.62668, 0.32331],[-0.92511, 0.32506 ], ],
    [ [ -0.61557, 0.32310],[-0.61582, 0.17344],[-0.31658, 0.17219],[-0.31641, 0.32200],[-0.61557, 0.32310 ], ],
    [ [ -0.30509, 0.32217],[-0.30522, 0.17238],[-0.00533, 0.17192],[-0.00528, 0.32179],[-0.30509, 0.32217 ], ],
    [ [ 0.00613, 0.32189],[0.00603, 0.17202],[0.30591, 0.17216],[0.30593, 0.32196],[0.00613, 0.32189 ], ],
    [ [ 0.31747, 0.32198],[0.31748, 0.17218],[0.61671, 0.17309],[0.61664, 0.32274],[0.31747, 0.32198 ], ],
    [ [ 0.62784, 0.32277],[0.62791, 0.17311],[0.92629, 0.17481],[0.92627, 0.32435],[0.62784, 0.32277 ], ],
    [ [ -1.08035, 0.16192],[-1.08040, 0.01245],[-0.78210, 0.00992],[-0.78203, 0.15954],[-1.08035, 0.16192 ], ],
    [ [ -0.77072, 0.15965],[-0.77083, 0.01002],[-0.47205, 0.00825],[-0.47194, 0.15807],[-0.77072, 0.15965 ], ],
    [ [ -0.46056, 0.15788],[-0.46069, 0.00805],[-0.16102, 0.00706],[-0.16092, 0.15700],[-0.46056, 0.15788 ], ],
    [ [ -0.14938, 0.15712],[-0.14941, 0.00717],[0.15061, 0.00715],[0.15061, 0.15711],[-0.14938, 0.15712 ], ],
    [ [ 0.16175, 0.15688],[0.16173, 0.00693],[0.46138, 0.00772],[0.46136, 0.15755],[0.16175, 0.15688 ], ],
    [ [ 0.47294, 0.15757],[0.47295, 0.00776],[0.77175, 0.00934],[0.77172, 0.15898],[0.47294, 0.15757 ], ],
    [ [ 0.78288, 0.15914],[0.78293, 0.00951],[1.08118, 0.01197],[1.08117, 0.16145],[0.78288, 0.15914 ], ],
    [ [ -1.08051, -0.00221],[-1.08062, -0.15168],[-0.78229, -0.15458],[-0.78222, -0.00496],[-1.08051, -0.00221 ], ],
    [ [ -0.77090, -0.00489],[-0.77092, -0.15452],[-0.47215, -0.15653],[-0.47212, -0.00671],[-0.77090, -0.00489 ], ],
    [ [ -0.46081, -0.00680],[-0.46076, -0.15662],[-0.16114, -0.15760],[-0.16114, -0.00767],[-0.46081, -0.00680 ], ],
    [ [ -0.14973, -0.00771],[-0.14978, -0.15765],[0.15021, -0.15775],[0.15028, -0.00781],[-0.14973, -0.00771 ], ],
    [ [ 0.16141, -0.00798],[0.16146, -0.15793],[0.46109, -0.15683],[0.46107, -0.00701],[0.16141, -0.00798 ], ],
    [ [ 0.47258, -0.00693],[0.47247, -0.15680],[0.77124, -0.15501],[0.77134, -0.00541],[0.47258, -0.00693 ], ],
    [ [ 0.78279, -0.00540],[0.78263, -0.15503],[1.08092, -0.15257],[1.08103, -0.00311],[0.78279, -0.00540 ], ],
    [ [ -0.92571, -0.16812],[-0.92560, -0.31760],[-0.62719, -0.32007],[-0.62733, -0.17043],[-0.92571, -0.16812 ], ],
    [ [ -0.61609, -0.17058],[-0.61594, -0.32023],[-0.31677, -0.32182],[-0.31685, -0.17203],[-0.61609, -0.17058 ], ],
    [ [ -0.30535, -0.17203],[-0.30526, -0.32183],[-0.00546, -0.32238],[-0.00545, -0.17251],[-0.30535, -0.17203 ], ],
    [ [ 0.00599, -0.17273],[0.00593, -0.32257],[0.30574, -0.32209],[0.30586, -0.17229],[0.00599, -0.17273 ], ],
    [ [ 0.31721, -0.17229],[0.31708, -0.32209],[0.61624, -0.32060],[0.61642, -0.17095],[0.31721, -0.17229 ], ],
    [ [ 0.62766, -0.17093],[0.62742, -0.32060],[0.92587, -0.31820],[0.92603, -0.16874],[0.62766, -0.17093 ], ],
    [ [ -0.92568, -0.33228],[-0.92551, -0.48171],[-0.62701, -0.48423],[-0.62725, -0.33471],[-0.92568, -0.33228 ], ],
    [ [ -0.61620, -0.33497],[-0.61596, -0.48448],[-0.31691, -0.48617],[-0.31702, -0.33654],[-0.61620, -0.33497 ], ],
    [ [ -0.30532, -0.33665],[-0.30529, -0.48627],[-0.00565, -0.48704],[-0.00552, -0.33739],[-0.30532, -0.33665 ], ],
    [ [ 0.00557, -0.33746],[0.00555, -0.48717],[0.30521, -0.48654],[0.30538, -0.33684],[0.00557, -0.33746 ], ],
    [ [ 0.31703, -0.33683],[0.31682, -0.48646],[0.61590, -0.48489],[0.61620, -0.33538],[0.31703, -0.33683 ], ],
    [ [ 0.62745, -0.33525],[0.62719, -0.48474],[0.92568, -0.48229],[0.92589, -0.33285],[0.62745, -0.33525 ], ],
    [ [ -0.77053, -0.49796],[-0.77039, -0.64731],[-0.47166, -0.64969],[-0.47181, -0.50032],[-0.77053, -0.49796 ], ],
    [ [ -0.46055, -0.50025],[-0.46036, -0.64962],[-0.16113, -0.65089],[-0.16118, -0.50150],[-0.46055, -0.50025 ], ],
    [ [ -0.14968, -0.50146],[-0.14964, -0.65087],[0.14989, -0.65100],[0.15006, -0.50159],[-0.14968, -0.50146 ], ],
    [ [ 0.16131, -0.50174],[0.16111, -0.65115],[0.46031, -0.65015],[0.46069, -0.50077],[0.16131, -0.50174 ], ],
    [ [ 0.47199, -0.50057],[0.47168, -0.64992],[0.77041, -0.64783],[0.77070, -0.49848],[0.47199, -0.50057 ], ],
    [ [ -0.61537, -0.66335],[-0.61517, -0.81261],[-0.31628, -0.81441],[-0.31649, -0.66519],[-0.61537, -0.66335 ], ],
    [ [ -0.30527, -0.66522],[-0.30512, -0.81439],[-0.00586, -0.81513],[-0.00583, -0.66592],[-0.30527, -0.66522 ], ],
    [ [ 0.00563, -0.66608],[0.00566, -0.81526],[0.30490, -0.81455],[0.30506, -0.66537],[0.00563, -0.66608 ], ],
    [ [ 0.31632, -0.66540],[0.31612, -0.81459],[0.61504, -0.81292],[0.61531, -0.66368],[0.31632, -0.66540 ], ],
    [ [ -0.46006, -0.82837],[-0.46002, -0.97759],[-0.16106, -0.97892],[-0.16098, -0.82981],[-0.46006, -0.82837 ], ],
    [ [ -0.14966, -0.82964],[-0.14964, -0.97869],[0.14946, -0.97881],[0.14963, -0.82976],[-0.14966, -0.82964 ], ],
    [ [ 0.16100, -0.82983],[0.16083, -0.97889],[0.45984, -0.97803],[0.46012, -0.82884],[0.16100, -0.82983 ], ],
];

decam_ccd_names = [ 'S29', 'S30', 'S31',
                        'S25', 'S26', 'S27', 'S28',
                        'S20', 'S21', 'S22', 'S23', 'S24',
                        'S14', 'S15', 'S16', 'S17', 'S18', 'S19',
                        'S8',  'S9',  'S10', 'S11', 'S12', 'S13',
                        'S1',  'S2',  'S3',  'S4',  'S5',  'S6',  'S7',
                        'N1',  'N2',  'N3',  'N4',  'N5',  'N6',  'N7',
                        'N8',  'N9',  'N10', 'N11', 'N12', 'N13',
                        'N14', 'N15', 'N16', 'N17', 'N18', 'N19',
                        'N20', 'N21', 'N22', 'N23', 'N24',
                        'N25', 'N26', 'N27', 'N28',
                        'N29', 'N30', 'N31',
                      ];

def draw_decam(ax,ra,dec,**kwargs):
    '''
    Helper function to draw the outline of DECam CCD's centered at a point
    
    Args:
      ax: matplotlib Axes; axes to draw CCD's over
      ra: float; Ra (in degrees) to center on
      dec: float; Dec (in degrees) to center on
      kwargs: dict; dictionary passed to ax.plot
    
    Returns:
      None
    
    '''
    decam_arr = np.array(decam_rds)
    for i in range(62):
        x_corner = decam_arr[i,:,0]
        y_corner = decam_arr[i,:,1]
        ax.plot(ra + x_corner,dec + y_corner/np.cos(np.pi * dec/180),transform = ax.get_transform('icrs'),**kwargs)

# run if not initialized from import (I'll be using these functions a lot, so this will be useful)
if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise Exception("Proper Usage: python draw_skymap.py REPO CLN")
    
    repo_name = sys.argv[1]
    cln = sys.argv[2]
    
    # query ned for pointing
    ned_result = Ned.query_object(cln)
    ra_cl = ned_result[0]['RA']
    dec_cl = ned_result[0]['DEC']
    
    butler = Butler(repo_name)

    # retrieve the skymap from the butler
    skymap = butler.get('skyMap',collections='skymaps',dataId={'instrument':'DECam','tract':0,'skymap':'{CLN}_skymap'.format(CLN=cln)})
    
    # a lazy trick to get the 'tract' object, which contains info about patches
    for tract in skymap:
        print(skymap.config)

    # creating the wcs
    wcs = WCS(tract.getWcs().getFitsMetadata())
    ax = pl.subplot(projection=wcs)
    ax.tick_params(axis='both',which='major',labelsize=8)
    
    # getting patch info and drawing
    id_array,bbox_array = get_patch_info(tract)
    draw_patches(ax,id_array,bbox_array,textparams={'horizontalalignment':'center','fontsize':'x-small'})
    
    # draw decam fov
    # in-case it's useful I'll leave this in... but the 1deg circle does a slightly better job
    # draw_decam(ax,ra_cl,dec_cl,linestyle='--',color='green',alpha=0.3,linewidth=0.5)
    
    # draw effective decam fov
    r = SphericalCircle((ra_cl*u.deg,dec_cl*u.deg), 1 * u.deg,
                     edgecolor='black', facecolor='none',
                     transform=ax.get_transform('icrs'),linestyle='--',
                     alpha=0.7,linewidth=0.5)
    ax.add_patch(r)
    
    # query hips for legacy-image
    hips_wcs, result = get_cluster_cutout(cln,fov=2)
    
    # transform and display the inverted image on the plot
    # flip axis to match fits wcs-convention
    ax.imshow(1-np.flip(result,axis=0),transform=ax.get_transform(hips_wcs))
    
    pl.savefig("skymap_indices.png",dpi=1080,bbox_inches='tight')
    
    
    
    





