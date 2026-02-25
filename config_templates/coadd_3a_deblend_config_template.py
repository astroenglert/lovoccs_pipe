config.saveLogOutput=True

# Minimum number of iterations before the optimizer is allowed to stop.
# default is 15
config.multibandDeblend.minIter=15

# Maximum number of iterations to deblend a single parent
# default is 300, let's crank this up to several thousand
config.multibandDeblend.maxIter=300

# Change in the loss function between iterations to exit fitter. Typically this is `1e-2` if measurements will be made on the flux re-distributed models and `1e-4` when making measurements on the models themselves.
config.multibandDeblend.relativeError=0.01

# Fraction of background RMS a pixel must haveto be included in the initial morphology
config.multibandDeblend.morphThresh=1.0

# The version of scarlet to use.
config.multibandDeblend.version='lite'

# The optimizer to use for fitting parameters and is only used when version='lite'
config.multibandDeblend.optimizer='adaprox'

# The type of image to use for initializing the morphology. Must be either 'chi2' or 'wavelet'. 
config.multibandDeblend.morphImage='chi2'

# Fraction of background to use for a sparsity threshold. This prevents sources from growing unrealistically outside the parent footprint while still modeling flux correctly for bright sources.
config.multibandDeblend.backgroundThresh=0.25

# Maximum number of proximal operator iterations inside of each iteration of the optimizer. This config field is only used if version='lite' and optimizer='adaprox'.
config.multibandDeblend.maxProxIter=1

# Number of wavelet scales to use for wavelet initialization. This field is only used when `version`='lite' and `morphImage`='wavelet'.
config.multibandDeblend.waveletScales=5

# Whether or not use use inverse variance weighting.If `useWeights` is `False` then flat weights are used
config.multibandDeblend.useWeights=True

# Model PSF side length in pixels
config.multibandDeblend.modelPsfSize=11

# Define sigma for the model frame PSF
config.multibandDeblend.modelPsfSigma=0.8

# Minimum Signal to noise to accept the source.Sources with lower flux will be initialized with the PSF but updated like an ordinary ExtendedSource (known in scarlet as a `CompactSource`).
config.multibandDeblend.minSNR=50.0

# Whether or not to save the SEDs and templates
config.multibandDeblend.saveTemplates=True

# Whether or not to process isolated sources in the deblender
config.multibandDeblend.processSingles=True

# Type of convolution to render the model to the observations.
# - 'fft': perform convolutions in Fourier space
# - 'real': peform convolutions in real space.
config.multibandDeblend.convolutionType='fft'

# Whether or not to solve for the best-fit spectra during initialization. This makes initialization slightly longer, as it requires a convolution to set the optimal spectra, but results in a much better initial log-likelihood and reduced total runtime, with convergence in fewer iterations.This option is only used when peaks*area < `maxSpectrumCutoff` will use the improved initialization.
config.multibandDeblend.setSpectra=True

# Whether or not to process isolated sources in the deblender
config.multibandDeblend.badMask=['BAD', 'CR', 'NO_DATA', 'SAT', 'SUSPECT', 'EDGE']

# Mask planes to ignore when performing statistics
config.multibandDeblend.statsMask=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended. If the fraction is `0` then the limit is a single pixel.
config.multibandDeblend.maskLimits={}

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
# default is 200, BCG can get close to enclosing this many sources...
config.multibandDeblend.maxNumberOfPeaks=200

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
# default is 100000... but that's not quite enough to get the bcg, which can extend ~5'x5' boxes
# let's be a bit more generous than that and try 3e6, for 1500px x 1500px
config.multibandDeblend.maxFootprintArea=int(3e6)

# Maximum rectangular footprint area * nPeaks in the footprint. This was introduced in DM-33690 to prevent fields that are crowded or have a LSB galaxy that causes memory intensive initialization in scarlet from dominating the overall runtime and/or causing the task to run out of memory. (<= 0: unlimited)

# default is 10000000, again this won't quite be enough. Assuming maximum footprint, let's let it deblend
# up to 200 sources
config.multibandDeblend.maxAreaTimesPeaks=int(6e8)

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
# default is zero, let's make sure it's disabled
config.multibandDeblend.maxFootprintSize=-1

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
# default is zero, let's make sure it's disabled
config.multibandDeblend.minFootprintAxisRatio=-1

# Maximum number of pixels * number of sources in a blend. This is different than `maxFootprintArea` because this isn't the footprint area but the area of the bounding box that contains the footprint, and is also multiplied by the number ofsources in the footprint. This prevents large skinny blends with a high density of sources from running out of memory. If `maxSpectrumCutoff == -1` then there is no cutoff.

# default is 1000000 this disables spectrum utilization... when obscured by the BCG how much benefit
# is there to that anyways? let's update it to 3e7 to keep it one order of mag above the maxarea
# like with the default settings
config.multibandDeblend.maxSpectrumCutoff=int(3e7)

# Whether or not to fallback to a smaller number of components if a source does not initialize
config.multibandDeblend.fallback=True

# Mask name for footprints not deblended, or None
config.multibandDeblend.notDeblendedMask='NOT_DEBLENDED'

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.multibandDeblend.catchFailures=True

# Columns to pass from the parent to the child. The key is the name of the column for the parent record, the value is the name of the column to use for the child.
config.multibandDeblend.columnInheritance={'deblend_nChild': 'deblend_parentNChild', 'deblend_nPeaks': 'deblend_parentNPeaks', 'deblend_spectrumInitFlag': 'deblend_spectrumInitFlag', 'deblend_blendConvergenceFailedFlag': 'deblend_blendConvergenceFailedFlag'}

# Names of flags which should never be deblended.
config.multibandDeblend.pseudoColumns=['merge_peak_sky', 'sky_source']

# Limit the number of sources deblended for CI to prevent long build times
config.multibandDeblend.useCiLimits=False

# Only deblend parent Footprints with a number of peaks in the (inclusive) range indicated.If `useCiLimits==False` then this parameter is ignored.
config.multibandDeblend.ciDeblendChildRange=[5, 10]

# Only use the first `ciNumParentsToDeblend` parent footprints with a total peak count within `ciDebledChildRange`. If `useCiLimits==False` then this parameter is ignored.
config.multibandDeblend.ciNumParentsToDeblend=10

# Identifier for a data release or other version to embed in generated IDs. Zero is reserved for IDs with no embedded release identifier.
config.idGenerator.release_id=0

# Number of (contiguous, starting from zero) `release_id` values to reserve space for. One (not zero) is used to reserve no space.
config.idGenerator.n_releases=1

# Mapping from band name to integer to use in the packed ID. The default (None) is to use a hard-coded list of common bands; pipelines that can enumerate the set of bands they are likely to see should override this.
config.idGenerator.packer.bands=None

# Number of bands to reserve space for. If zero, bands are not included in the packed integer at all. If `None`, the size of 'bands' is used.
config.idGenerator.packer.n_bands=0

# Number of tracts, or, more precisely, one greater than the maximum tract ID.Default (None) obtains this value from the skymap dimension record.
config.idGenerator.packer.n_tracts=None

# Number of patches per tract, or, more precisely, one greater than the maximum patch ID.Default (None) obtains this value from the skymap dimension record.
config.idGenerator.packer.n_patches=None

# name for connection inputSchema
config.connections.inputSchema='{inputCoaddName}Coadd_mergeDet_schema'

# name for connection peakSchema
config.connections.peakSchema='{inputCoaddName}Coadd_peak_schema'

# name for connection mergedDetections
config.connections.mergedDetections='{inputCoaddName}Coadd_mergeDet'

# name for connection coadds
config.connections.coadds='{inputCoaddName}Coadd_calexp'

# name for connection outputSchema
config.connections.outputSchema='{outputCoaddName}Coadd_deblendedFlux_schema'

# name for connection fluxCatalogs
config.connections.fluxCatalogs='{outputCoaddName}Coadd_deblendedFlux'

# name for connection templateCatalogs
config.connections.templateCatalogs='{outputCoaddName}Coadd_deblendedModel'

# name for connection deblendedCatalog
config.connections.deblendedCatalog='{outputCoaddName}Coadd_deblendedCatalog'

# name for connection scarletModelData
config.connections.scarletModelData='{outputCoaddName}Coadd_scarletModelData'

# Template parameter used to format corresponding field template parameter
config.connections.inputCoaddName='deep'

# Template parameter used to format corresponding field template parameter
config.connections.outputCoaddName='deep'

