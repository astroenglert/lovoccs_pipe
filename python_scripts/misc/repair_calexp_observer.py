from lsst.daf.butler import Butler
import astropy.io.fits as fits
import sys

# iterate through dataset_type in collection and find all headers containing a single apostraphe
def repair_headers(butler,dataset_type='calexp',collection='DECam/processing/calexp_r'):
    '''
    Iterate through a dataset_type and find/replace all headers containing an apostraphe (which crashes LSP at later steps)
    
    Args:
      butler: Butler; butler for the repository
      dataset_type: string; datasetType to load, defaults to 'calexp'
      collection: string; collection to load datasets from, defaylts to 'DECam/processing/calexp_r'
      
    Returns:
      None
    
    '''
    # write the output to disk?
    write = False
    
    # run the query to get calexp datasetrefs
    query = butler.registry.queryDatasets(dataset_type,collections=collection)
    
    # iterate through query results
    for file in query:
    
        # get the filepaths and hdul
        filepath = butler.getURI(file)
        hdul = fits.open(filepath.ospath)
        
        # check for, and replace, any apostraphe's (the problem is actually repeated apostraphe's I think
        # but better safe than sorry)
        obs = hdul[0].header['OBSERVER']
        prp = hdul[0].header['PROPOSER']
        dtp = hdul[0].header['DTPI']
        hdul.close()
        
        if "'" in obs:
            print(f'Updating observer for {filepath}')
            fits.setval(filepath.path,keyword='OBSERVER',value=obs.replace("'",""))
        if "'" in prp:
            print(f'Updating proposer for {filepath}')
            fits.setval(filepath.path,keyword='PROPOSER',value=prp.replace("'",""))
        if "'" in dtp:
            print(f'Updating drpi for {filepath}')
            fits.setval(filepath.path,keyword='DRPI',value=dtp.replace("'",""))
        
    return

# if ran from cln, execute the following
if __name__ == "__main__":
    
    # checking for correct usage
    if len(sys.argv)!=5:
        print("python this_script.py cln repo dataset_type collection")
        sys.exit(1)
    
    # collect input arguments
    cln = sys.argv[1]
    butler = Butler(sys.argv[2])
    dataset_type = sys.argv[3]
    collection = sys.argv[4]
    
    # update the headers
    repair_headers(butler,dataset_type,collection)
    
    
    




