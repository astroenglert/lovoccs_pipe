# query NED to get the cluster coordinates
from astroquery.ipac.ned import Ned

ned_result = Ned.query_object('cluster_name')
ra_cl = ned_result[0]['RA']
dec_cl = ned_result[0]['DEC']
zL = ned_result[0]['Redshift']

config.skyMap["discrete"].projection = "TAN" # tangent-projection
config.skyMap["discrete"].pixelScale = 0.263 # DECam angular resolution

# inner-patch size, these tile the tract exactly
config.skyMap["discrete"].patchInnerDimensions = [4000,4000] 

# outer border extends +100px beyond inner-patch (but are truncated at tract borders), creates a net 200px overlap region 
config.skyMap["discrete"].patchBorder = 100

# I'm going to FIX the number of patches and center it on the cluster

config.skyMap["discrete"].raList = [ra_cl]
config.skyMap["discrete"].decList = [dec_cl]
config.skyMap["discrete"].radiusList = [0.7]
config.skyMap.name = "discrete"

config.name = "cluster_name_skymap"

