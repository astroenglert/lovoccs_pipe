# Estimate the background again after final source detection?
config.detection.reEstimateBackground=True

# type of statistic to use for grid points
config.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
# default is 4096, trying to decrease it to 128x128
config.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
# and let's measure that background with the good ol' Chebyshev method
config.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detection.background.weighting=True

