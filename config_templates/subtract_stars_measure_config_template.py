### Configuration for task `measureExtendedPsf'

# Size, in pixels, of the subregions over which the stacking will be iteratively performed.
config.stack_bright_stars.subregion_size=[100, 100]

# Apply magnitude cut before stacking?
config.stack_bright_stars.do_mag_cut=False

# Magnitude limit, in Gaia G; all stars brighter than this value will be stacked
config.stack_bright_stars.mag_limit=12.0

# Mapping from detector IDs to focal plane region names. If empty, a constant extended PSF model is built from all selected bright stars.
# detectors need multiple stars to build an accurate model... so I'm grouping them a little bit

#TODO fix these to be detector:"region", v27 is "region":[det_ids] but I'm on v26 rn
#config.detectors_focal_plane_regions={
#    1:[1,2,3],
#    2:[4,8,9],
#    3:[5,6,10],
#    4:[7,11,12],
#    5:[13,14,19,20],
#    6:[15,16,21,22],
#    7:[17,18,23,24],
#    8:[25,26,32,33],
#    9:[27,28,29,34,35,36],
#    10:[30,31,37,38],
#    11:[39,40,45,46],
#    12:[41,42,47,48],
#    13:[43,44,49,50],
#    14:[51,52,56],
#    15:[53,57,58],
#    16:[54,55,59],
#    17:[60,61,62],
#}

config.detectors_focal_plane_regions={
    1:'z',2:'z',3:'z',
    4:'o1',8:'o1',9:'o1',
    5:'o2',6:'o2',10:'o2',
    7:'o3',11:'o3',12:'o3',
    13:'tw1',14:'tw1',19:'tw1',20:'tw1',
    15:'tw2',16:'tw2',21:'tw2',22:'tw2',
    17:'tw3',18:'tw3',23:'tw3',24:'tw3',
    25:'th1',26:'th1',32:'th1',33:'th1',
    27:'th2',29:'th2',29:'th2',34:'th2',35:'th2',36:'th2',
    30:'th3',31:'th3',37:'th3',38:'th3',
    39:'fo1',40:'fo1',45:'fo1',46:'fo1',
    41:'fo2',42:'fo2',47:'fo2',48:'fo2',
    43:'fo3',44:'fo3',49:'fo3',50:'fo3',
    51:'fi1',52:'fi1',56:'fi1',
    53:'fi2',57:'fi2',58:'fi2',
    54:'fi3',55:'fi3',59:'fi3',
    60:'s',61:'s',62:'s',
}

# name for connection input_brightStarStamps
config.connections.input_brightStarStamps='brightStarStamps'

# name for connection extended_psf
config.connections.extended_psf='extended_psf'
