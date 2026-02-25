### Configuration for task `processBrightStars'
# Flag to enable/disable saving of log output for a task, enabled by default.
config.saveLogOutput=True

# Magnitude limit, in Gaia G; all stars brighter than this value will be processed.
config.magLimit=12.0

# Size of the stamps to be extracted, in pixels.
config.stampSize=[1000, 1000]

# 'Buffer' factor to be applied to determine the size of the stamp the processed stars will be saved in. This will also be the size of the extended PSF model.
config.modelStampBuffer=1.1

# Apply transform to bright star stamps to correct for optical distortions?
config.doApplyTransform=True

# Inner and outer radii of the annulus used to compute AnnularFlux for normalization, in pixels.
#TODO Find the best annuli for this
config.annularFluxRadii=[180, 220]

# Minumum number of valid pixels that must fall within the annulus for the bright star to be saved for subsequent generation of a PSF.
config.minValidAnnulusFraction=0.0

# Apply full focal plane sky correction before extracting stars?
config.doApplySkyCorr=True

# Should stars with NaN annular flux be discarded?
config.discardNanFluxStars=False

# Padding to add to 4 all edges of the bounding box (pixels)
config.refObjLoader.pixelMargin=250

# Always use this reference catalog filter, no matter whether or what filter name is supplied to the loader. Effectively a trivial filterMap: map all filter names to this filter. This can be set for purely-astrometric catalogs (e.g. Gaia DR2) where there is only one reasonable choice for every camera filter->refcat mapping, but not for refcats used for photometry, which need a filterMap and/or colorterms/transmission corrections.
config.refObjLoader.anyFilterMapsToThis=None

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist in the refcat. Note that this does not perform any bandpass corrections: it is just a lookup.
config.refObjLoader.filterMap={}

# Require that the fields needed to correct proper motion (epoch, pm_ra and pm_dec) are present?
config.refObjLoader.requireProperMotion=False

# name for connection skyCorr
config.connections.skyCorr='skyCorr'

# name for connection refCat
# hard-coded default
config.connections.refCat='gaia'

# name for connection brightStarStamps
config.connections.brightStarStamps='brightStarStamps'
