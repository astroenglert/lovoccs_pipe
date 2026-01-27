# a small python script for storing filters to be used by mass_map
#TODO do we have to change these to operate on-sky rather than pixels (flat-sky approximation)?

import numpy as np


def schirmer_filter(radius,aperture_size=8000,x_cut=0.15,a=6,b=150,c=47,d=50,*_):
    '''
    The Schirmer Filter, a filter which is optimized for detecting NFW-like structures in shear-fields.
    
    Args:
        radius: Numpy array; an array of radii to evaluate the filter on
        aperture_size: float-like; the 'schirmer-radius' of the filter
        x_cut: float-like; specifies the filter-sloap and sets the characteristic-scale of the filter to x_cut*smoothing
    
    Returns:
        Q; Numpy array; an array containing the filter evaluated at each radius
    
    '''
    
    x = radius/aperture_size
    Q = ( 1/( 1 + np.exp(a - b*x) + np.exp(-c + d*x)) )*( np.tanh(x/x_cut)/(x/x_cut) )
    Q = Q/(np.pi * aperture_size**2)
    
    return Q


def truncated_gaussian(radius,aperture_size=8000,*_):
    '''
    Some teams prefer using a truncated-gaussian filter (e.g. Miyazaki+14, Hamana+12, among many others).
    
    Args:
        radius: Numpy array; an array of radii to evaluate the filter on
        smoothing: float-like; the 'radius' of the filter
    
    Returns:
        Q; Numpy array; an array containing the filter evaluated at each radius
    
    '''
    
    evaluate_gaussian = radius < aperture_size
    truncated = ~evaluate_gaussian
    
    gauss = lambda x : 1 - (1 + x**2)*np.exp(-x**2)
    
    x = radius/aperture_size
    Q = np.piecewise(x, [evaluate_gaussian,truncated], [gauss,lambda x: 0] )
    Q = Q/(np.pi * radius**2)
    
    return Q


def polynomial_filter(radius,aperture_size=8000,l=1,*_):
    '''
    A family of polynomial filters, not-commonly used but easy to implement (Schneider+96)
    
    Args:
        radius: Numpy array; an array of radii to evaluate the filter on
        smoothing: float-like; the 'radius' of the filter
        l: float-like; a constant specifying the filter
        
    Returns:
        Q; Numpy array; an array containing the filter evaluated at each radius
    
    '''
    
    evaluate_polynomial = radius < aperture_size
    truncated = ~evaluate_polynomial
    
    poly = lambda x : (x**2) * (1 - x**2)**2
    
    x = radius/aperture_size
    Q = np.piecewise(x, [evaluate_polynomial,truncated], [poly,lambda x: 0] )
    Q = (l+1)*(l+2)*Q/(np.pi * aperture_size**2)
    
    return Q
    

#TODO this one is a bit of a pain bc it requires solving for a/b/c based on continuity
# nice project for an undergrad(?) to tinker-with and see how it works relative to the above
def truncated_sis(radius,aperture_size=8000,nu1=0.1,nu2=0.3,*_):
    '''
    A filter optimized for detecting SIS, but truncated, is also common (Oguri+21, Schneider+96)
    
    Args:
        radius: Numpy array; an array of radii to evaluate the filter on
        smoothing: float-likwe; the 'radius' of the filter
        nu1: float-like; one of two constants which specify the filter
        nu2: float-like; one of two constants which specify the filter
        
    Returns:
        Q; Numpy array; an array containing the filter evaluated at each radius
    
    '''
    
    raise NotImplementedError("TODO")
    #return Q


# a dictionary of the implemented filters
implemented_filters = {'schirmer':schirmer_filter,'gaussian':truncated_gaussian,'polynomial':polynomial_filter}

