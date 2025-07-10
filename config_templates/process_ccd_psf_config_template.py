import lsst.meas.extensions.piff.piffPsfDeterminer

config.psfIterations=psf_iter
config.measurePsf.psfDeterminer.name='psf_determiner'

config.measurePsf.starSelector['objectSize'].fluxMin=min_flux
config.measurePsf.starSelector['objectSize'].widthStdAllowed=std_width
config.measurePsf.starSelector['objectSize'].doFluxLimit=flux_lim

# let the flux-cut take care of the low-SN stars
config.measurePsf.starSelector['objectSize'].doSignalToNoiseLimit=False
