from abc import ABC
import math

from importlib import resources as impresources
from pathlib import Path

import numpy as np

from scipy.special import gamma
from scipy.interpolate import interp1d

from astropy.table import Table
from astropy.io import ascii

# homebrew modules below
from . import SED_CWWSB, SED_COSMOS, SED_TEST

sed_choice = SED_COSMOS
sed_db = impresources.files(sed_choice)

class TemplateBase(ABC):

    # initialize the template
    def __init__(self,filter_transmissions,redshifts = np.arange(0.01,1.5,0.001) ,wavelengths = np.linspace(1e-3,3e4,int(3e4)),cosmology=None,):
        
        # setting the cosmology accordingly
        # optional, is here if we want to use surface-brightness or experiment with other more complex templates/priors
        self.cosmo = cosmology
        
        # setting the redshifts to evaluate along
        self.redshifts = redshifts
        
        # setting the wavelengths to evaluate along
        self.wavelengths = wavelengths
        
        # loading the transmission
        self.filter_transmissions = filter_transmissions
        
        # load the sed from disk
        self.sed_data = np.loadtxt(self.sed_filepath)
        
        # build an interpolated sed
        self.sed = self._interp_sed(self.sed_data)
        
        # there are proper caching libraries that we need to test, but defining variables is good enough for now
        self._flux_cache = None
        
        pass
    
    # return a function for the prior, e.g. p(T|m0)
    # write this to return either a float or an array of len(m0), numpy-friendly! 
    def get_template_prior(self,m0):
        raise NotImplementedError
    
    # get prior for the redshift, e.g. p(z|T,m0)
    # write this to return either an array of len(z) or an array of (len(z),len(m0)), numpy-friendly!
    def get_redshift_prior(self,m0,redshift=None):
        raise NotImplementedError
    
    # the location of the SED on-disk
    def sed_filepath(self):
        raise NotImplementedError
    
    # name of the SED, just in-case
    def sed_name(self):
        raise NotImplementedError
    
    # run a simple linear-interpolation on the SED
    def _interp_sed(self,sed_data):
        
        interpolated_sed = interp1d(x=sed_data[:,0],y=sed_data[:,1],bounds_error=False,fill_value=0)
        
        return interpolated_sed    
    
    # convenience function for hubble-like priors
    def _hubble_redshift_prior_amplitude(self,alpha_t, z_0t, k_mt, m0, m_min):
        tmp = (1./(k_mt*(m0-m_min)+z_0t) )**(1.+alpha_t) *alpha_t**2 / gamma(1./alpha_t)
        return tmp
    
    # redshift-priors derived for Hubble for three spectral classes (Coe+06)
    def _hubble_redshift_priors(self,m0,galaxy_type,z=None):
        
        if z is None:
            z = self.redshifts
        
        if not (galaxy_type in [1,2,3]):
            print(f'No galaxy of type {galaxy_type}')
            raise Exception('No galaxy of the provided type')
        
        # this is all carried-over straight from Shenming's old photo-z script
        m_min = 20
        m0_new = np.where(m0 < m_min, m_min, m0)
        
        if galaxy_type == 1:
            alpha_t = 2.465
            z_0t = 0.431
            k_mt = 0.0913
            Amp = self._hubble_redshift_prior_amplitude(alpha_t, z_0t, k_mt, m0_new, m_min)
            result = Amp[:,None] * (z**alpha_t+0.) * np.exp( -( z / ( z_0t + k_mt*(m0_new[:,None] - m_min) ))**alpha_t )
            
            # experimental gaussian prior at mean-redshift of the cluster
            # currently disabled
            sigma = 0.03
            weight = 0.5
            cluster_z = 0.093
            tmp = np.exp(-(z-cluster_z)**2/(2.*sigma**2))/(sigma*(2.*np.pi)**0.5)
            result = result * weight + tmp * (1.-weight)
        
        elif galaxy_type == 2:
            alpha_t = 1.806
            z_0t = 0.390
            k_mt = 0.0636
            Amp = self._hubble_redshift_prior_amplitude(alpha_t, z_0t, k_mt, m0_new, m_min)
            result = Amp[:,None] * (z**alpha_t+0.) * np.exp( -( z / ( z_0t + k_mt*(m0_new[:,None]-m_min) ))**alpha_t )
        
        else:
            alpha_t = 0.906
            z_0t = 0.0626
            k_mt = 0.123
            Amp = self._hubble_redshift_prior_amplitude(alpha_t, z_0t, k_mt, m0_new, m_min)
            result = Amp[:,None] * z**alpha_t * np.exp( -( z / ( z_0t + k_mt*(m0_new[:,None]-m_min) ))**alpha_t )

        return result
        
    # template-priors derived for Hubble for three spectral classes (Coe+06)
    def _hubble_template_priors(self,m0,galaxy_type):
    
        if not (galaxy_type in [1,2,3]):
            print(f'No galaxy of type {galaxy_type}')
            raise Exception('No galaxy of the provided type')
        
        m_min = 20
        delta_m = np.where(m0 < m_min, 0, m0 - m_min)
        
        early_ft = 0.35
        early_kt = 0.45
        early_types = early_ft * np.exp( -early_kt * delta_m )
        
        late_ft = 0.50
        late_kt = 0.147
        late_types = late_ft * np.exp( -late_kt * delta_m )
        
        other_types = 1 - late_types - early_types
        
        if galaxy_type == 1:
            return early_types
        elif galaxy_type == 2:
            return late_types
        else:
            return other_types
    
    # priors developed by folks running the Next Generation Virgo Cluster Survey
    # reaches similar depths (m~25, SN~10) and reduces outliers at m<20 (low-z)
    def _ngvcs_template_priors(self,m0,galaxy_type):
        return np.ones(len(m0))
        
        if not (galaxy_type in [1,2,3]):
            print(f'No galaxy of type {galaxy_type}')
            raise Exception('No galaxy of the provided type')
        
        # use np.piecewise to break up m0 into the appropriate magnitude ranges
        condition_1 = m0 < 12.5
        condition_2 = (m0 >= 12.5) & (m0 < 17)
        condition_3 = (m0 >= 17) & (m0 < 20)
        condition_4 = (m0 >= 20)
        conditions = [condition_1,condition_2,condition_3,condition_4]
        
        ref_mag = [0,12.5,17,20]
        
        early_ft = [1/3,0.86,0.65,0.30]
        early_kt = [0,0.062,0.257,0.40]
        
        late_ft = [1/3,0.14,0.23,0.35]
        late_kt = [0,-0.108,-0.014,0.30]
        
        # lazy trick, use arrays of lambda-functions for storing these
        early_prior = [1] * 4
        late_prior = [1] * 4
        for i in range(4):
            early_prior[i] = lambda x : early_ft[i] * np.exp( - early_kt[i] * ( x - ref_mag[i] ) )
            late_prior[i] = lambda x : late_ft[i] * np.exp( - late_kt[i] * ( x - ref_mag[i] ) )
        
        prior = np.zeros(len(m0))
        if galaxy_type == 1:
            for i in range(4):
                prior[conditions[i]] = early_prior[i](m0[conditions[i]])

                # # experimental gaussian fit
                # # add_gaussian(z, prior) returns a new prior?
                # sigma = 0.03
                # weight = 0.5 # fixed or vary by magnitude 
                # cluster_z = 0.093
                # tmp = np.exp(-(z-cluster_z)**2/(2.*sigma**2))/(sigma*(2.*np.pi)**0.5)
                # prior = prior * weight + tmp * (1.-weight)
            #prior = np.piecewise(m0,conditions,early_prior)
            return prior
        elif galaxy_type == 2:
            for i in range(4):
                prior[conditions[i]] = late_prior[i](m0[conditions[i]])
            return prior
        else:
            for i in range(4):
                prior[conditions[i]] = 1 - late_prior[i](m0[conditions[i]]) - early_prior[i](m0[conditions[i]])
            return prior
        
    def _ngvcs_redshift_priors(self,m0,galaxy_type,z=None):
        
        if z is None:
            z = self.redshifts
        
        return np.ones((len(m0),len(z)))

        if not (galaxy_type in [1,2,3]):
            print(f'No galaxy of type {galaxy_type}')
            raise Exception('No galaxy of the provided type')
        
        # use np.piecewise to break up m0 into the appropriate magnitude ranges
        condition_1 = m0 < 12.5
        condition_2 = (m0 >= 12.5) & (m0 < 17)
        condition_3 = (m0 >= 17) & (m0 < 20)
        condition_4 = (m0 >= 20)
        conditions = [condition_1,condition_2,condition_3,condition_4]
        
        ref_mag = [0,12.5,17,20]
        
        early_alpha = [1,2.46,2.46,2.46]
        early_z0 = [1,0,0.121,0.431]
        early_kmt = [1,0.027,0.103,0.091]
        
        late_alpha = [1,2.07,1.94,1.81]
        late_z0 = [1,0,0.095,0.390]
        late_kmt = [1,0.021,0.098,0.100]
        
        irr_alpha = [1,1.89,1.95,2.00]
        irr_z0 = [1,0,0.069,0.300]
        irr_kmt = [1,0.15,0.077,0.150]
        
        early_prior = [1] * 4
        early_amp = [1] * 4
        late_prior = [1] * 4
        late_amp = [1] * 4
        irr_prior = [1] * 4
        irr_amp = [1] * 4
        for i in range(4):
            early_prior[i] = lambda m0 : ( z**(early_alpha[i]) ) * np.exp( - ( z/(early_z0[i] + early_kmt[i]*(m0[:,None] - ref_mag[i]) ) )**(early_alpha[i]) )
            late_prior[i] = lambda m0 : ( z**(late_alpha[i]) ) * np.exp( - ( z/(late_z0[i] + late_kmt[i]*(m0[:,None] - ref_mag[i]) ) )**(late_alpha[i]) )
            irr_prior[i] = lambda m0 : ( z**(irr_alpha[i]) ) * np.exp( - ( z/(irr_z0[i] + irr_kmt[i]*(m0[:,None] - ref_mag[i]) ) )**(irr_alpha[i]) )
        
        # manually override the first prior to be flat
        # this isn't doing what it;s supposed too, but ends up as a flat constant in the end so it doesn't actually matter (by a happy coincidence!)
        early_prior[0] = lambda x : np.ones((len(x),len(z)))/(np.max(z) - np.min(z))
        late_prior[0] = lambda x : np.ones((len(x),len(z)))/(np.max(z) - np.min(z))
        irr_prior[0] = lambda x : np.ones((len(x),len(z)))/(np.max(z) - np.min(z))
        
        # can't use np.interpolate since the outputs are 2D arrays, manually use the conditions
        output_arr = np.zeros((len(m0),len(z)))
        
        if galaxy_type == 1:
            for i in range(4):
                amp = self._hubble_redshift_prior_amplitude(early_alpha[i], early_z0[i], early_kmt[i], m0[conditions[i]], ref_mag[i])
                output_arr[conditions[i]] = amp[:,None]*early_prior[i](m0[conditions[i]])
            
        elif galaxy_type == 2:
            for i in range(4):
                amp = self._hubble_redshift_prior_amplitude(late_alpha[i], late_z0[i], late_kmt[i], m0[conditions[i]], ref_mag[i])
                output_arr[conditions[i]] = late_prior[i](m0[conditions[i]])
                
        else:
            for i in range(4):
                amp = self._hubble_redshift_prior_amplitude(irr_alpha[i], irr_z0[i], irr_kmt[i], m0[conditions[i]], ref_mag[i])
                output_arr[conditions[i]] = irr_prior[i](m0[conditions[i]])
                
        return output_arr
    
    def _cosmos_template_priors(self,m0,galaxy_type):
        output_arr = np.ones(len(m0))
        return output_arr

    def _cosmos_redshift_priors(self,m0,galaxy_type,z=None):
        if z is None:
            z = self.redshifts
        output_arr = np.ones((len(m0),len(z)))
        return output_arr

    # get the flux along the redshifts for this template
    # how should this interface with the cache? Need to decide
    #TODO somewhere there is a bug (some missing factor or something) in the code computing fluxes here that I introduced after computing the cached fluxes; this needs to be patched before trying to generate new cached fluxes
    def compute_flux(self,filter_transmissions=None,wavelength=None,redshift=None):
        
        if redshift is None:
            redshift = self.redshifts
        if wavelength is None:
            wavelength = self.wavelengths
        if filter_transmissions is None:
            filter_transmissions = self.filter_transmissions
        
        fluxes = Table()
        for filter,trans in self.filter_transmissions.transmission_functions.items():
            
            result = np.trapz(
                              y = self.sed(wavelength[:,None]/(1 + redshift))/(1.+redshift)*trans(wavelength[:,None])*wavelength[:,None],
                              x = wavelength,
                              axis=0
                             )
            result = result / self.filter_transmissions.get_normalization(filter,wavelength)
            fluxes[filter] = result
            fluxes['redshift'] = redshift
        
        # store these in the "cache"
        self._flux_cache = fluxes
        
        return fluxes
    
    # write computed redshifts/fluxes to the disk
    def write_flux_to_disk(self,filepath,**kwargs):
        
        if self._flux_cache is None:
            print("Flux cache is empty... now computing fluxes")
            fluxes = self.compute_flux(**kwargs)
        
        ascii.write(self._flux_cache,filepath,format='csv',overwrite=True)
        return True
    
    # draw the sed
    def draw_template(self,ax,wavelengths=None,**kwargs):
        
        if wavelengths is None:
            wavelengths = self.wavelengths
        
        ax.plot(wavelengths,self.sed(wavelengths),label=self.sed_name,**kwargs)
        return ax

    # load redshifts/fluxes from disk
    def load_flux_from_disk(self,filepath,filter_transmission):
        
        fluxes = ascii.read(filepath,format='csv')
        self._flux_cache = fluxes
        self.redshifts = fluxes['redshift']
        
        return True

def make_template_class(sed_name: str, template_type: int):
    """Ning 9/16 Factory for making TemplateBase subclasses."""
    print(sed_name, template_type)
    # sed_path = sed_db.joinpath(f"{sed_name}.sed")
    _sed_name = sed_name

    class TemplateSubclass(TemplateBase):
        sed_filepath = sed_db.joinpath(f"{_sed_name}.sed")
        sed_name = _sed_name

        def get_template_prior(self, m0):
            return self._cosmos_template_priors(m0, template_type)

        def get_redshift_prior(self, m0, redshift=None):
            return self._cosmos_redshift_priors(m0, template_type, redshift)

    TemplateSubclass.__name__ = sed_name 
    return TemplateSubclass

sed_list = []
with impresources.open_text(sed_choice, "cosmossedswdust136.list") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name = Path(line).stem  # strip directory + extension
        sed_list.append(name)

cosmos_sed_list = [make_template_class(name, idx) for idx, name in enumerate(sed_list, start=1)]

class ElB2004a(TemplateBase):
    
    sed_filepath = sed_db.joinpath('El_B2004a.sed')
    #sed_filepath = sed_db.joinpath('LePhare_Templates/CWW_E_ext.sed')
    
    sed_name = 'EL_B2004a'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,1)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,1,redshift)
    
class SbcB2004a(TemplateBase):
    
    sed_filepath = sed_db.joinpath('Sbc_B2004a.sed')
    #sed_filepath = sed_db.joinpath('LePhare_Templates/CWW_Sbc_ext.sed')
    
    sed_name = 'Sbc_B2004a'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,2)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,2,redshift)

class ScdB2004a(TemplateBase):
    
    sed_filepath = sed_db.joinpath('Scd_B2004a.sed')
    #sed_filepath = sed_db.joinpath('LePhare_Templates/CWW_Scd_ext.sed')
    
    sed_name = 'Scd_B2004a'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,2)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,2,redshift)

class ImB2004a(TemplateBase):
    
    sed_filepath = sed_db.joinpath('Im_B2004a.sed')
    #sed_filepath = sed_db.joinpath('LePhare_Templates/CWW_Im_ext.sed')
    
    sed_name = 'Im_B2004a'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,3)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,3,redshift)

class SB2B2004a(TemplateBase):
    
    sed_filepath = sed_db.joinpath('SB2_B2004a.sed')
    #sed_filepath = sed_db.joinpath('LePhare_Templates/KIN_SB2_ext.sed')
    
    sed_name = 'SB2_B2004a'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,3)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,3,redshift)

class SB3B2004a(TemplateBase):
    
    sed_filepath = sed_db.joinpath('SB3_B2004a.sed')
    #sed_filepath = sed_db.joinpath('LePhare_Templates/KIN_SB3_ext.sed')
    
    sed_name = 'SB3_B2004a'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,3)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,3,redshift)

class ssp25MyrZ008(TemplateBase):
    
    sed_filepath = sed_db.joinpath('ssp_25Myr_z008.sed')
    
    sed_name = 'ssp_25Myr_z008'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,3)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,3,redshift)

class ssp5MyrZ008(TemplateBase):
    
    sed_filepath = sed_db.joinpath('ssp_5Myr_z008.sed')
    
    sed_name = 'ssp_5Myr_z008'
    
    def get_template_prior(self,m0):
        return self._ngvcs_template_priors(m0,3)
    
    def get_redshift_prior(self,m0,redshift=None):
        return self._ngvcs_redshift_priors(m0,3,redshift)

# convenience variable containing our default list of templates
# ordered based on galaxy-evolution
default_sed_list = [
                    ElB2004a,
                    SbcB2004a,
                    ScdB2004a,
                    ImB2004a,
                    SB2B2004a,
                    SB3B2004a,
                    ssp25MyrZ008,
                    ssp5MyrZ008,
                   ]

# dictionary of sed-list for given prior sets
prior_dict = {'ngvs':default_sed_list,'cosmos':cosmos_sed_list}

