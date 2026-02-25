### Configuration for task `starSub'
# Magnitude limit, in Gaia G; all stars brighter than this value will be subtracted
config.magLimit=12.0

# How the model should be scaled to each bright star; implemented options are `annularFlux` to reuse the annular flux of each stamp, or `leastSquare` to perform least square fitting on each pixel with no bad mask plane set.
config.scalingType='leastSquare'

# Apply full focal plane sky correction before extracting stars?
config.doApplySkyCorr=True

# name for connection inputExposure
config.connections.inputExposure='calexp'

# name for connection inputBrightStarStamps
config.connections.inputBrightStarStamps='brightStarStamps'

# name for connection inputExtendedPsf
config.connections.inputExtendedPsf='extended_psf'

# name for connection skyCorr
config.connections.skyCorr='skyCorr'

# name for connection outputExposure
config.connections.outputExposure='{outputExposureName}_calexp'

# name for connection outputBackgroundExposure
config.connections.outputBackgroundExposure='{outputBackgroundName}_calexp_background'

# Template parameter used to format corresponding field template parameter
config.connections.outputExposureName='brightStar_subtracted'

# Template parameter used to format corresponding field template parameter
config.connections.outputBackgroundName='brightStars'

