# currently not implemented, but may be worth experimenting with as an alt. to check/select-visit

# Lower bound of quantile range to select. Sorts visits by seeing from narrow to wide, and select those in the interquantile range (qMin, qMax). Set qMin to 0 for Best Seeing. This config should be changed from zero only for exploratory diffIm testing.
config.qMin=0.0

# Upper bound of quantile range to select. Sorts visits by seeing from narrow to wide, and select those in the interquantile range (qMin, qMax). Set qMax to 1 for Worst Seeing.
config.qMax=0.33

# At least this number of visits selected and supercedes quantile. For example, if 10 visits cover this patch, qMin=0.33, and nVisitsMin=5, the best 5 visits will be selected.
config.nVisitsMin=6

# Do remove visits that do not actually overlap the patch?
config.doConfirmOverlap=True

# Minimum visit MJD to select
config.minMJD=None

# Maximum visit MJD to select
config.maxMJD=None

# name for connection skyMap
config.connections.skyMap='skyMap'

# name for connection visitSummaries
config.connections.visitSummaries='finalVisitSummary'

# name for connection goodVisits
config.connections.goodVisits='goodSeeingVisits'

# Template parameter used to format corresponding field template parameter
config.connections.coaddName='goodSeeing'

# Template parameter used to format corresponding field template parameter
config.connections.calexpType=''
