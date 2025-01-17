import sys

import numpy as np

from astropy.table import Table

#TODO we should do something more sophisticated eventually, SED-fitting?
def separate_stars_galaxies(catalog,extendedness_flag='extendedness'):
    '''
    Separate stars from galaxies in catalog according to the extendedness_flag. I'm defining this as a separate function since there are other algorithms and ways of doing this that may be worth exploring (outside of simply cmodel-psf magnitude, which is the default)
    
    Args:
        catalog: Astropy Table; table containing extendedness_flag to separate galaxies
        extendedness_flag: String; defaults to 'extendedness'
        
    Returns:
        catalog_stars: Astropy Table; copy of catalog containing just stars
        catalog_gals: Astropy Table; copy of catalog containing just galaxies
    '''
    
    is_star = catalog['extendedness'] < 0.5
    is_gal = catalog['extendedness'] > 0.5
    
    catalog_stars = catalog[is_star]
    catalog_gals = catalog[is_gal]
    
    return catalog_stars, catalog_gals


if __name__ == '__main__':
    
    if len(sys.argv)!=4:
        print("Usage: python separate_stars_galaxies.py catalog_all catalog_star catalog_galaxy")
        sys.exit(1)

    catalog_all_filename = sys.argv[1]
    catalog_star_filename = sys.argv[2]
    catalog_gal_filename = sys.argv[3]
    
    # load the catalog
    catalog_all = Table.read(catalog_all_filename, format="ascii.csv")
    catalog_stars,catalog_gals = separate_stars_galaxies(catalog_all,extendedness_flag='extendedness')
    
    # write the star/gal catalogs
    catalog_stars.write(catalog_star_filename, format="ascii.csv", overwrite=True)
    catalog_gals.write(catalog_gal_filename, format="ascii.csv", overwrite=True)

