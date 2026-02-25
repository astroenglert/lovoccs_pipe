# need to override the default reference catalog connection

# given that there are ~50visits/band in the best cases (lots of archival data), 25th quantile is reasonable to fill in the gaps btwn ccd's with dithering.
# default is 33rd
config.qMax=0.25

# 10 isn't quite enough to cover up oversaturated regions or ccd-gaps, let's try 15 instead
config.nVisitsMin=15

# might be useful in the future to have these available
# defaults are None
config.minMJD=None
config.maxMJD=None
