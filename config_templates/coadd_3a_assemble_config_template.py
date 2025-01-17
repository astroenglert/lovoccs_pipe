# options we may need/want to tweak in the future

# in case we want to try psfMatched
config.warpType='direct'

# mean v median stacking
config.statistic='MEAN'

# in-case we need to tweak clipping
config.sigmaClip=3.0

# select-visit for good-seeing coadd
config.doSelectVisits=False

# we've already selected the best visits, so I'm disabling any further cuts
config.select.maxEllipResidual=None
config.select.maxScaledSizeScatter=None
config.select.maxPsfTraceRadiusDelta=None

# need to disable these for the static-sky model (artifact rejection)
config.assembleStaticSkyModel.select.maxEllipResidual=None
config.assembleStaticSkyModel.select.maxScaledSizeScatter=None
config.assembleStaticSkyModel.select.maxPsfTraceRadiusDelta=None

# minimum to pick out transients without significantly masking persistent sources at patch border
config.maxNumEpochs=3

# debugging obj detection
# causes a type-mismatch somewhere which causes the pipeline to file, so no debug logs for me :'(
# config.assembleStaticSkyModel.doWrite=True

# this flag is a little different in gen2 and the 'suspect' pixel region may be correlated with my artifact, DEFAULT: ['NO_DATA', 'BAD', 'SAT', 'SUSPECT']
# Looks like this config option is fine, the real issue was the psf-cuts above
# config.badMaskPlanes=['NO_DATA', 'BAD', 'SAT', 'EDGE']

# these will be useful for debugging
config.inputRecorder.saveErrorCcds=True
config.inputRecorder.saveEmptyCcds=True
