# a small python script for storing profiles to fit with mass_fit
# each of these functions returns the convergence and shears (gamma1/gamma2)

#TODO we should re-write this to be object oriented, start with a 'radially_symmetric_profile' base-class
# then the following classes define the methods to compute mean_sigma(r) and sigma(r)

import numpy as np

# default fiducial cosmology
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


# helper function to be used in mass_fit, not called explicitly in methods below to save a bit of runtime
def sigma_crit(zS,zL=0.05,cosmo=cosmo):
    '''
    Computes the critical surface-density at a given zS,zL, and cosmology.
    
    Args:
        zS: Numpy array; an (N,)-array containing the redshift of background galaxies
        zL: float; a float specifying the redshift of the cluster
        cosmo: Astropy Cosmology; a reference cosmology
        
    
    Returns:
        sigma_crit: Numpy array; an (N,)-array containing the critical surface density for a given background galaxy
        
    '''
    
    distance_ratio = cosmo.angular_diameter_distance(zS) / (cosmo.angular_diameter_distance(zL)*cosmo.angular_diameter_distance_z1z2(zL,zS))
    sigma_crit = ( ((3e5)**2)/(4*np.pi*4.3011790220362e-09) ) * distance_ratio.value
    
    return sigma_crit


# convergance/shears for an nfw-profile
# analytic forms derived from Wright+00
# mass-concentration relation from Child+18
#TODO generalize to different mass-concentration relations
def shear_nfw(M200,delta_x,delta_y,zS=None,zL=0.05,cosmo=cosmo):
    '''
    An NFW profile, the current-standard for modelling the density profile of halos.
    
    Args:
        M200: float; the mass (in Msol) at a radius enclosing 200 x mean-overdensity
        zS: Numpy array; an (N,)-array containing the redshift of background galaxies. If None this returns delta_sigma
        delta_x: Numpy array; an (N,)-array containing the x-dist (in radians) between the background galaxy and cluster-center
        delta_y: Numpy array; an (N,)-array containing the y-dist (in radians) between the background galaxy and cluster-center
        zL: float; a float specifying the redshift of the cluster
        cosmo: Astropy Cosmology; a reference cosmology
    
    
    Returns:
        kappa; Numpy array; an (N,)-array containing the convergance at each background-galaxy
        gamma_1; Numpy array; an (N,)-array containing the 'x-comp' of shear at each background-galaxy
        gamma_2; Numpy array; an (N,)-array containing the 'y-comp' of shear at each background-galaxy
    
    '''
    
    # first compute the concentration and the critical overdensity
    concentration = (75.4) * ( (1+zL)**(-0.422) ) * (M200**(-0.089))
    delta_crit = (200/3) * (concentration**3)/(np.log(1+concentration) - (concentration/(1+concentration)))
    
    # next get the critical-density in units of Msol/Mpc^3
    rho_crit = cosmo.critical_density(zL).to(u.solMass/(u.Mpc**3)).value
    
    # now we can compute r200 and the scale-radius, units of Mpc
    r_200 = ( (3/(4*np.pi)) * (M200 / (200*rho_crit)))**(1/3)
    r_scale = r_200/concentration
    
    # last of the pre-requisites, compute the critical surface density
    if zS is not None:
        distance_ratio = cosmo.angular_diameter_distance(zS) / (cosmo.angular_diameter_distance(zL)*cosmo.angular_diameter_distance_z1z2(zL,zS))
        sigma_crit = ( ((3e5)**2)/(4*np.pi*4.3011790220362e-09) ) * distance_ratio.value
    
    # define the dimensionless parameter-x, project into the lens-plane by scaling by DL
    x = np.sqrt(delta_x**2 + delta_y**2)*cosmo.angular_diameter_distance(zL).value / r_scale
        
    #TODO we should reformat this code, use np.interpolate or something to manage the different cases
    xx = np.abs(x)
    Sigma = np.zeros_like(xx)
    g_tmp = np.zeros_like(xx)
    
    # idx: index
    # 0<X<1
    idxa = (xx>0.0)
    idxb = (xx<1.0)
    idx1 = idxa&idxb
    x1 = 1.0/(xx[idx1]*xx[idx1]-1.0)
    x2_up = 2.0*np.arctanh(np.sqrt((1.0-xx[idx1])/(1.0+xx[idx1])))
    x2_down = np.sqrt(1.0-xx[idx1]*xx[idx1])
    x2 = x2_up/x2_down
    x3 = np.log(xx[idx1]/2.0)
    x4 = 4.0/(xx[idx1]*xx[idx1])
    Sigma[idx1] = 2.0*r_scale*delta_crit*rho_crit*x1*(1.0-x2)
    g_tmp[idx1] = x2*x4 + x4*x3 - 2.0*x1 + 2.0*x2*x1
    
    # X=1
    idx2 = (xx==1.0)
    Sigma[idx2] = 2.0*r_scale*delta_crit*rho_crit/3.0
    g_tmp[idx2] = 10.0/3.0 + 4.0*np.log(0.5)

    # X>1
    idx3 = (xx>1.0)
    x1 = 1.0/(xx[idx3]*xx[idx3]-1.0)
    x2_up = 2.0*np.arctan(np.sqrt((xx[idx3]-1.0)/(1.0+xx[idx3])))
    x2_down = np.sqrt(xx[idx3]*xx[idx3]-1.0)
    x2 = x2_up/x2_down
    x3 = np.log(xx[idx3]/2.0)
    x4 = 4.0/(xx[idx3]*xx[idx3])
    Sigma[idx3] = 2.0*r_scale*delta_crit*rho_crit*x1*(1.0-x2)
    g_tmp[idx3] = x2*x4 + x4*x3 - 2.0*x1 + 2.0*x2*x1
    
    # if the source redshift isn't specified, spit-out delta_sigma
    if zS is None:
        return g_tmp * r_scale * delta_crit * rho_crit
    
    kappa = Sigma/sigma_crit

    gamma = r_scale*delta_crit*rho_crit*g_tmp/sigma_crit

    
    # now the shears can be computed easily
    theta = np.arctan2(delta_y,delta_x)
    gamma_1 = - gamma*np.cos(2*theta)
    gamma_2 = - gamma*np.sin(2*theta)
    
    return kappa, gamma_1, gamma_2

#TODO we should implement this functionality for an SIS at a later-date
def shear_sis(M200,zS,delta_x,delta_y,zL=0.05,cosmo=cosmo):
    '''
    An SIS profile, falls by 1/r^2.
    
    Args:
        M200: float; the mass (in Msol) at a radius enclosing 200 x mean-overdensity
        zS: Numpy array; an (N,)-array containing the redshift of background galaxies. If None, return delta_Sigma.
        delta_x: Numpy array; an (N,)-array containing the x-dist (in arcsec) between the background galaxy and cluster-center
        delta_y: Numpy array; an (N,)-array containing the y-dist (in arcsec) between the background galaxy and cluster-center
        zL: float; a float specifying the redshift of the cluster
        cosmo: Astropy Cosmology; a reference cosmology
        
    
    Returns:
        kappa; Numpy array; an (N,)-array containing the convergance at each background-galaxy
        gamma_1; Numpy array; an (N,)-array containing the 'x-comp' of shear at each background-galaxy
        gamma_2; Numpy array; an (N,)-array containing the 'y-comp' of shear at each background-galaxy
    
    '''
    
    raise NotImplementedError


#TODO we should implement this functionality for at a later date
def shear_elliptical_nfw(M200,zS,delta_x,delta_y,zL=0.05,cosmo=cosmo,ecc=1):
    '''
    An elliptical nfw profile.
    
    Args:
        M200: float; the mass (in Msol) at a radius enclosing 200 x mean-overdensity
        ecc: float; the eccentricity parameter of the profile
        zS: Numpy array; an (N,)-array containing the redshift of background galaxies
        delta_x: Numpy array; an (N,)-array containing the x-dist (in arcsec) between the background galaxy and cluster-center
        delta_y: Numpy array; an (N,)-array containing the y-dist (in arcsec) between the background galaxy and cluster-center
        zL: float; a float specifying the redshift of the cluster
        cosmo: Astropy Cosmology; a reference cosmology
        
    
    Returns:
        kappa; Numpy array; an (N,)-array containing the convergance at each background-galaxy
        gamma_1; Numpy array; an (N,)-array containing the 'x-comp' of shear at each background-galaxy
        gamma_2; Numpy array; an (N,)-array containing the 'y-comp' of shear at each background-galaxy
    
    '''
    
    raise NotImplementedError



#TODO we should implement this functionality for at a later date
def shear_elliptical_sis(M200,zS,delta_x,delta_y,zL=0.05,cosmo=cosmo,ecc=1):
    '''
    An elliptical SIS profile.
    
    Args:
        M200: float; the mass (in Msol) at a radius enclosing 200 x mean-overdensity
        ecc: float; the eccentricity parameter of the profile
        zS: Numpy array; an (N,)-array containing the redshift of background galaxies
        delta_x: Numpy array; an (N,)-array containing the x-dist (in arcsec) between the background galaxy and cluster-center
        delta_y: Numpy array; an (N,)-array containing the y-dist (in arcsec) between the background galaxy and cluster-center
        zL: float; a float specifying the redshift of the cluster
        cosmo: Astropy Cosmology; a reference cosmology
        
    
    Returns:
        kappa; Numpy array; an (N,)-array containing the convergance at each background-galaxy
        gamma_1; Numpy array; an (N,)-array containing the 'x-comp' of shear at each background-galaxy
        gamma_2; Numpy array; an (N,)-array containing the 'y-comp' of shear at each background-galaxy
    
    '''
    
    raise NotImplementedError


# defining a list of the implemented shear_profiles
implemented_profiles = {'nfw':shear_nfw}

    
