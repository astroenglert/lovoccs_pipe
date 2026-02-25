# need to override the default reference catalog connection

config.connections.refCat="photom_ref"

config.match.refObjLoader.filterMap={ 'process_band' : 'photom_mag' }

# may need to update the name of this dataset type
# config.connections.sourceTableHandles='sourceTable_visit'
