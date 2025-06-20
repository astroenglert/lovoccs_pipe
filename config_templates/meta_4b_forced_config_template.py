
# running a stripped-down version of measurments for metadetection
config.measurement.plugins.names=['base_CircularApertureFlux', 'base_GaussianFlux', 'base_InputCount', 'base_LocalBackground', 'base_PixelFlags', 'base_PsfFlux', 'base_SdssCentroid', 'base_SdssShape', 'base_TransformedCentroid', 'base_TransformedShape', 'base_Variance', 'modelfit_CModel', 'modelfit_DoubleShapeletPsfApprox']

# only run the essential aperture flux
config.measurement.plugins['base_CircularApertureFlux'].radii=[12.0]

# run undeblended measurements necessary
config.measurement.undeblended.names=['base_PsfFlux']

# allow apcorr for psf_flux
config.doApCorr=True
config.applyApCorr.doFlagApCorrFailures=True

# allow the essential aperture corrections
config.applyApCorr.proxies={'undeblended_base_PsfFlux_instFlux': 'base_PsfFlux_instFlux'}

