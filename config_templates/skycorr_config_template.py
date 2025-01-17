# Alright, the algorithm here is:
#
# 1. Detect and mask sources
# 2. Subtract a large-scale background (bgmodel1) to remove gradients
# 3. Subtract a sky-frame to remove large-scale fixed patterns (sky)
# 4. Subtract a small-scale background (bgmodel2) to remove small-scale residuals
# 5. Iterate 1-4 3x (default), in theory this preserves LSB features with the right settings

# -- MASKING SETTINGS -- #

# Number of masking iterations
config.maskObjects.nIter=3

# These options control source masking prior to background subtraction
# the "subtractBackground" options here don't actually create a background used later... but is
# an intermediate step in the masking process 
 
# how large a region of the sky should be used for each background point
config.maskObjects.subtractBackground.binSize=1024

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.subtractBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.subtractBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.subtractBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.subtractBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.maskObjects.subtractBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.maskObjects.subtractBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.maskObjects.subtractBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.maskObjects.subtractBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.maskObjects.subtractBackground.weighting=True

# detected sources with fewer than the specified number of pixels will be ignored
config.maskObjects.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.maskObjects.detection.isotropicGrow=False

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.maskObjects.detection.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.maskObjects.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.maskObjects.detection.returnOriginalFootprints=False

# We'll want to be very generous with the detection-threshold of footprints
# Fewer sigma's (default is 2.5) create more generous masks, this option
# might need a lot of tweaking

# Threshold for footprints; exact meaning and units depend on thresholdType.
config.maskObjects.detection.thresholdValue=2.5

# Include threshold relative to thresholdValue
config.maskObjects.detection.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
config.maskObjects.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
config.maskObjects.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.maskObjects.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.maskObjects.detection.reEstimateBackground=False

# type of statistic to use for grid points
config.maskObjects.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.maskObjects.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.maskObjects.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.maskObjects.detection.background.weighting=True

# type of statistic to use for grid points
config.maskObjects.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.detection.tempLocalBackground.binSizeY=0

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.maskObjects.detection.nPeaksMaxSimple=1

# Mask planes to exclude when detecting sources.
config.maskObjects.detection.excludeMaskPlanes=[]

# Detection threshold (standard deviations)
config.maskObjects.detectSigma=5.0

# Interpolate when removing objects?
config.maskObjects.doInterpolate=True

# type of statistic to use for grid points
config.maskObjects.interpolate.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.interpolate.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.interpolate.binSize=256

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.interpolate.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.interpolate.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.interpolate.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.interpolate.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Iteratively mask objects to find good sky?
config.doMaskObjects=True

# -- LARGE SCALE CONFIGS -- #

# the binning for the large-scale scale (which handles gradients) has been set to match
# the defaults used to produce the skyframe, covering two detectors (4096px * 0.015mm/px)
# This is technically a "full focal-plane" background which extends beyond single detectors
# hence being measured in mm

# Bin size in x
config.bgModel1.xSize=61.44

# Bin size in y
config.bgModel1.ySize=61.44

# Pixel size in same units as xSize/ySize
config.bgModel1.pixelSize=0.015

# Minimum fraction of bin size for good measurement
config.bgModel1.minFrac=0.1

# Mask planes to treat as bad
config.bgModel1.mask=['BAD', 'SAT', 'INTRP', 'DETECTED', 'DETECTED_NEGATIVE', 'EDGE', 'NO_DATA']

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.bgModel1.interpolation='AKIMA_SPLINE'

# Do smoothing?
config.bgModel1.doSmooth=False

# Smoothing scale, as a multiple of the bin size
config.bgModel1.smoothScale=2.0

# the large-scale backgrond is sampled from the focal-plane in 64x64 superpixels onto a ccd
# e.g. the large-scale background is created from the 61.44mm (4096px) superpixels, but when sampled
# back to a ccd is done in 64px bins (if my understanding is correct... it might be wrong `\_(=>)_/`)

# Binning to use for CCD background model (pixels)
config.bgModel1.binning=64

# this controls the scaling of the sky-frame relative to the image (e.g. img depth determines
# how bright the sky-frame should be, so it is upscaled/downscaled as needed)
# it's usually carried out by randomly sampling the sky and the exposure with clipping to pick
# average-background values, then are upscaled/downscaled to match

# k-sigma rejection iterations for sky scale
config.sky.skyIter=3

# k-sigma rejection threshold for sky scale
config.sky.skyRej=3.0

# type of statistic to use for grid points
config.sky.background.statistic='MEANCLIP'

# these bins just control super-pixels used for measuring the background level
# not for carrying out the correction.

# Superpixel size in x
config.sky.background.xBinSize=32

# Superpixel size in y
config.sky.background.yBinSize=32

# How to interpolate the background values. This maps to an enum; see afw::math::Background
config.sky.background.algorithm='NATURAL_SPLINE'

# Names of mask planes to ignore while estimating the background
config.sky.background.mask=['SAT', 'BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# -- SKYFRAME CONFIGS -- #

# Number of samples in x for scaling sky frame
config.sky.xNumSamples=4

# Number of samples in y for scaling sky frame
config.sky.yNumSamples=4

# type of statistic to use for grid points
config.sky.stats.statistic='MEANCLIP'

# Clipping threshold for background
config.sky.stats.clip=3.0

# Clipping iterations for background
config.sky.stats.nIter=3

# Mask planes to reject
config.sky.stats.mask=['SAT', 'DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'NO_DATA']

# Do sky frame subtraction?
config.doSky=True

# -- SMALL SCALE CONFIGS -- #

# after subtracting the skyframe... there is occassionally some smaller structures left behind
# so a second smaller-scale background model is fit. This is also a "focal plane model", so we
# have to enter things in mm

# defaults here are 256px (* 0.015mm/px), to catch background at a scale below the skyframe
# Bin size in x
config.bgModel2.xSize=3.84

# Bin size in y
config.bgModel2.ySize=3.84

# Pixel size in same units as xSize/ySize
config.bgModel2.pixelSize=0.015

# Minimum fraction of bin size for good measurement
config.bgModel2.minFrac=0.5

# Mask planes to treat as bad
config.bgModel2.mask=['BAD', 'SAT', 'INTRP', 'DETECTED', 'DETECTED_NEGATIVE', 'EDGE', 'NO_DATA']

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.bgModel2.interpolation='AKIMA_SPLINE'

# Do smoothing?
config.bgModel2.doSmooth=True

# Smoothing scale, as a multiple of the bin size
config.bgModel2.smoothScale=1.0

# Binning to use for CCD background model (pixels)
config.bgModel2.binning=64

# Do final (cleanup) background model subtraction, after sky frame subtraction?
config.doBgModel2=True

# -- THE REST CONTROL INPUT/OUTPUT CONNECTIONS -- #

# Binning factor for constructing full focal plane '*_camera' output datasets
config.binning=8

# name for connection rawLinker
config.connections.rawLinker='raw'

# name for connection calExps
config.connections.calExps='calexp'

# name for connection calBkgs
config.connections.calBkgs='calexpBackground'

# name for connection skyFrames
config.connections.skyFrames='sky'

# name for connection camera
config.connections.camera='camera'

# name for connection skyCorr
config.connections.skyCorr='skyCorr'

# name for connection calExpMosaic
config.connections.calExpMosaic='calexp_skyCorr_visit_mosaic'

# name for connection calBkgMosaic
config.connections.calBkgMosaic='calexpBackground_skyCorr_visit_mosaic'

