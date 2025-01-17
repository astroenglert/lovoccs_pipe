# config for finalCharacterization

# we might want to tweak some of these parameters at a later stage...

# by default this order is 2
config.psf_determiner['piff'].spatialOrder=piff_order

# default here is 25
config.psf_determiner['piff'].stampSize=piff_stamp

# save debug data from piff in-case we need to check piff
config.psf_determiner['piff'].debugStarData=True

# default SN for aperture correction is quite strict at 200, maybe tweak it
config.measure_ap_corr.sourceSelector['science'].signalToNoise.minimum=ap_corr_sn
