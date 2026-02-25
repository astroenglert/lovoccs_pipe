# this is a goodSeeing stack, default is 'deep'
# I think these are connections to the previous deepCoadd warps/psf-matched exps
# in which case these should be unchanged!
#config.coaddName='goodSeeing'
#config.select.connections.coaddName='goodSeeing'

# our internal cuts should take care of these already!
config.select.maxEllipResidual=None
config.select.maxScaledSizeScatter=None
config.select.maxPsfTraceRadiusDelta=None

# need to disable these for the static-sky model (artifact rejection)
config.assembleStaticSkyModel.select.maxEllipResidual=None
config.assembleStaticSkyModel.select.maxScaledSizeScatter=None
config.assembleStaticSkyModel.select.maxPsfTraceRadiusDelta=None

# minimum to pick out transients without significantly masking persistent sources at patch border
config.maxNumEpochs=3

# these will be useful for debugging
config.inputRecorder.saveErrorCcds=True
config.inputRecorder.saveEmptyCcds=True
