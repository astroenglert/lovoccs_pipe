
import lsst.meas.extensions.simpleShape

# disable this, no need to run any matching between coadds and calexps
config.doMatchSources=False
config.doPropagateFlags=False
config.doWriteMatchesDenormalized=False

# stripped down for measurement list for metadetect, only shape+photometry and other essetials
# HSM needs to be run to avoid a validation error, but no point in running others (e.g. KSB)
config.measurement.plugins.names=['base_Blendedness', 'base_CircularApertureFlux', 'base_GaussianFlux', 'base_InputCount', 'base_LocalBackground', 'base_LocalPhotoCalib', 'base_LocalWcs', 'base_PixelFlags', 'base_PsfFlux', 'base_SdssCentroid', 'base_SdssShape', 'base_SkyCoord', 'base_Variance', 'ext_shapeHSM_HsmPsfMoments', 'ext_shapeHSM_HsmShapeRegauss', 'ext_shapeHSM_HsmSourceMoments', 'ext_shapeHSM_HsmSourceMomentsRound', 'modelfit_CModel', 'ext_simpleShape_SimpleShape', 'modelfit_DoubleShapeletPsfApprox' ]

# 12.0 is the default used as a reference instFlux, which is required
config.measurement.plugins['base_CircularApertureFlux'].radii=[12.0]

# Maximum radius for pixels to include, in units of sigma
config.measurement.plugins['ext_simpleShape_SimpleShape'].nSigmaRegion=3.0

# Sigma of circular Gaussian used as weight function, in pixels
# chosen to match the FWHM of the larger PSF from metadetect, 1.2"
# chosen to match Sheldon+2023 (Metadetection for LSST)
config.measurement.plugins['ext_simpleShape_SimpleShape'].sigma=2.0

# Skip the aperture corrections
config.doApCorr=False

