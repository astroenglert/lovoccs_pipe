# using astrom_ref for astrometry
config.connections.astrometryRefCat = "astrom_ref"

# use a single band from the astrometry source catalog
config.astrometryRefObjLoader.anyFilterMapsToThis="astrom_band"

# some reference objects are missing pm-info
config.astrometryRefObjLoader.requireProperMotion=False

# using photom_ref for photometry
config.connections.photometryRefCat = "photom_ref"

# pointing LSP to the correct columns for each filter
# photom_ref has process_band magnitudes stored with header photom_mag
config.photometryRefObjLoader.filterMap = { 'process_band' : 'photom_mag' }

# we do not apply color-terms by default
config.applyColorTerms=False

# using polynomial of order jointcal_order for joint-calibration
config.photometryVisitOrder=jointcal_order

# it may be beneficial to adjust the astrometry order, these are the defaults for initial
config.astrometrySimpleOrder=3
config.astrometryChipOrder=1
config.astrometryVisitOrder=5

# by default, save the photom/astrom residuals in-case we need to debug
config.writeChi2FilesInitialFinal=True

# this directory has to exist beforehand, so create it during the setup step
config.debugOutputPath='jointcal_residuals/order_jointcal_order/'

# below are remnants from process_ccd_config_template.py, might need to enable these if jointcal crashes

# disable mag-lim since our catalogs already have this
# config.astrometry.referenceSelector.doMagLimit=False

# disable color-lims for the same reason
# config.photoCal.match.referenceSelection.colorLimits['g-r'].minimum=None
# config.photoCal.match.referenceSelection.colorLimits['g-r'].maximum=None
# config.photoCal.match.referenceSelection.colorLimits['r-i'].minimum=None
# config.photoCal.match.referenceSelection.colorLimits['r-i'].minimum=None

# Everything below has been hardcoded to point to r, g, etc. without using the filter map >:(

# despite mag & colorLims being disabled, these have to be specified to avoid a crash
# config.photoCal.match.referenceSelection.magLimit.fluxField='photom_mag_flux'
# config.photoCal.match.referenceSelection.colorLimits['g-r'].primary='photom_mag_flux'
# config.photoCal.match.referenceSelection.colorLimits['g-r'].secondary='photom_mag_flux'
# config.photoCal.match.referenceSelection.colorLimits['r-i'].primary='photom_mag_flux'
# config.photoCal.match.referenceSelection.colorLimits['r-i'].secondary='photom_mag_flux'

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
