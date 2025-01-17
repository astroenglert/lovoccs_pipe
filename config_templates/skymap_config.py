config.skyMap.projection = "TAN" # tangent-projection
config.skyMap.pixelScale = 0.263 # DECam angular resolution

# inner-patch size, these tile the tract exactly
config.skyMap.patchInnerDimensions = [4000,4000] 

# outer border extends +100px beyond inner-patch (but are truncated at tract borders), creates a net 200px overlap region 
config.skyMap.patchBorder = 100  

