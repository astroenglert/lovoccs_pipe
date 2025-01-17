# this is the replacement for makeCoaddTempExp
# currently not-implemented, but these seem to be the important configs we may want to tweak

# we produce a separate skycorr-coadded stack for LSB
config.doApplySkyCorr=True
config.bgSubtracted=True

# if we use star subtraction at some point this will need to be un-commented
#config.connections.calexpType='brightStar_subtracted_'

# for source injection this config option is going to be important
config.hasFakes=False

# we've already selected the best visits, so I'm disabling any future cuts
config.select.maxEllipResidual=None
config.select.maxScaledSizeScatter=None
config.select.maxPsfTraceRadiusDelta=None

# might help with debugging
config.inputRecorder.saveErrorCcds=True
config.inputRecorder.saveEmptyCcds=True
