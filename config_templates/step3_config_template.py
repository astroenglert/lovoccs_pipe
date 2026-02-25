# For now, I'm just going to record the config options we may need to tweak
# I'll note the task that is relevant and split the configs into separate files later

# FOR TESTING, start by running two visits to be co-added through the LSP

# == makeWarp == 
# this is the replacement for makeCoaddTempExp

# if we want to investigate the good-seeing coadd...
config.coaddName='deep'
config.select.connections.coaddName='deep'

# eventually we'll want to adjust this for Zach's stuff
config.doApplySkyCorr=False
config.bgSubtracted=True

# for source injection this config option is going to be important
config.hasFakes=False

# we should use direct coaddition (inverse-variance weighting)
# so I think we can disable this option... but need to test it
config.makePsfMatched=True

# == assembleCoadd == #

# in case we want to try psfMatched
config.warpType='direct'

# mean v median stacking
config.statistic='MEAN'

# in-case we need to tweak clipping
config.sigmaClip=3.0

# select-visit for good-seeing coadd
config.doSelectVisits=False

# == measure == #

config.coaddName='deep'

config.connections.refCat="photom_ref"

config.match.refObjLoader.filterMap={ 'process_band' : 'photom_mag' }

# == selectGoodSeeingVisits ==




# using astrom_ref for astrometry
config.connections.astromRefCat = "astrom_ref"

# use a single band from the astrometry source catalog
config.astromRefObjLoader.anyFilterMapsToThis="astrom_band"

# using photom_ref for photometry
config.connections.photoRefCat = "photom_ref"

# pointing LSP to the correct columns for each filter
# photom_ref has process_band magnitudes stored with header photom_mag
config.photoRefObjLoader.filterMap = { 'process_band' : 'photom_mag' }

# we do not apply color-terms by default
config.photoCal.applyColorTerms=False

# disable mag-lim since our catalogs already have this
config.astrometry.referenceSelector.doMagLimit=False

# disable color-lims for the same reason
config.photoCal.match.referenceSelection.colorLimits['g-r'].minimum=None
config.photoCal.match.referenceSelection.colorLimits['g-r'].maximum=None
config.photoCal.match.referenceSelection.colorLimits['r-i'].minimum=None
config.photoCal.match.referenceSelection.colorLimits['r-i'].minimum=None

# Everything below has been hardcoded to point to r, g, etc. without using the filter map >:(

# despite mag & colorLims being disabled, these have to be specified to avoid a crash
config.photoCal.match.referenceSelection.magLimit.fluxField='photom_mag_flux'
config.photoCal.match.referenceSelection.colorLimits['g-r'].primary='photom_mag_flux'
config.photoCal.match.referenceSelection.colorLimits['g-r'].secondary='photom_mag_flux'
config.photoCal.match.referenceSelection.colorLimits['r-i'].primary='photom_mag_flux'
config.photoCal.match.referenceSelection.colorLimits['r-i'].secondary='photom_mag_flux'

# these do not need to be specified if color-terms is disabled
# but, for clusters with complete ps1 coverage, we can apply the correction
# these configs are formatted, by default, for ps1

#config.photoCal.colorterms.data['ps1*'].data['g DECam SDSS c0001 4720.0 1520.0'].primary='gmag'
#config.photoCal.colorterms.data['ps1*'].data['g DECam SDSS c0001 4720.0 1520.0'].secondary='rmag'

#config.photoCal.colorterms.data['ps1*'].data['r DECam SDSS c0002 6415.0 1480.0'].primary='rmag'
#config.photoCal.colorterms.data['ps1*'].data['r DECam SDSS c0002 6415.0 1480.0'].secondary='imag'

#config.photoCal.colorterms.data['ps1*'].data['i DECam SDSS c0003 7835.0 1470.0'].primary='imag'
#config.photoCal.colorterms.data['ps1*'].data['i DECam SDSS c0003 7835.0 1470.0'].secondary='zmag'

#config.photoCal.colorterms.data['ps1*'].data['z DECam SDSS c0004 9260.0 1520.0'].primary='zmag'
#config.photoCal.colorterms.data['ps1*'].data['z DECam SDSS c0004 9260.0 1520.0'].secondary='ymag'
