from abc import ABC

from importlib import resources as impresources
from pathlib import Path

import numpy as np

from scipy.interpolate import interp1d

# homebrew modules below
from . import transmission_database

trans_db = impresources.files(transmission_database)

# rename this, maybe to FilterSetBase, nor sure
# by writing an abstract class, it'll be easier to expand all of our processing to different sets of filters
class TransmissionBase(ABC):
    
    def __init__(self):
    
        # initialize a list of filters
        self.filters = list(self.transmission_dictionary.keys())
        
        # read the transmission dictionary
        self.transmissions = self._load_transmission()
        
        # initialize interpolated transmission functions
        self.transmission_functions = self._interpolate_transmissions()
        
        pass
    
    # dictionary which stores filepaths to each filter transmission
    # two-columns, one wavelength in angstroms and another is the response
    @property
    def transmission_dictionary(self):
        raise NotImplementedError
    
    # string just for storing the instrument name
    @property
    def instrument_name(self):
        raise NotImplementedError

    # dictionary for storing mapping from filters to plotted colors
    @property
    def color_dict(self):
        raise NotImplementedError
       
    # loads transmissions according to transmission_dictionary
    def _load_transmission(self):
        
        transmission_dictionary = {}
        for filter, filename in self.transmission_dictionary.items():
            transmission_dictionary[filter] = np.loadtxt(filename)
        
        return transmission_dictionary
    
    # defines interpolated functions of the transmissions
    #TODO default is a simple linear interpolation, does a better inteprolation impact anything noticably?
    def _interpolate_transmissions(self):
        
        interpolated_transmissions = {}
        for filter, trans in self.transmissions.items():
            interpolated_transmissions[filter] = interp1d(x=trans[:,0],y=trans[:,1],bounds_error=False,fill_value=0)
            
        return interpolated_transmissions
    
    # computes the flux normalization-constant for a given band
    #TODO default is trapezoidal, should we do a bit better than this?
    def get_normalization(self,filter,wavelengths):
        
        # check that the filter is part of this set
        if not (filter in self.filters):
            print('f{filter} is not included in this class!')
            raise Exception("Filter not implemented!")
        
        # otherwise integrate over the range of wavelengths
        # int( T(Lambda) * c dLambda/Lambda )
        integrand = 1e6 * self.transmission_functions[filter](wavelengths) / wavelengths
        result = np.trapz(y=integrand,x=wavelengths)
        return result
    
    def draw_filter_transmissions(self,wavelengths,ax,**kwargs):
        
        for filter,trans in self.transmission_functions.items():
            ax.plot(wavelengths,trans(wavelengths),color=self.color_dict[filter],label=f'{filter}',**kwargs)
        
        return ax
        
# a class for the transmission of a few different instruments
# really this just neatly wraps-up constants for each set of filters
class DECamTransmission(TransmissionBase):
    
    transmission_dictionary = {
                               'u':trans_db.joinpath('u_DECam.res'),
                               'g':trans_db.joinpath('g_DECam.res'),
                               'r':trans_db.joinpath('r_DECam.res'),
                               'i':trans_db.joinpath('i_DECam.res'),
                               'z':trans_db.joinpath('z_DECam.res'),
                               'Y':trans_db.joinpath('Y_DECam.res'),
                              }
                              
    instrument_name = 'DECam'
    
    color_dict = {
                  'u':'C0',
                  'g':'C2',
                  'r':'C3',
                  'i':'C7',
                  'z':'C5',
                  'Y':'C6',
                 }
    
    pass

class HSCTransmission(TransmissionBase):
    
    transmission_dictionary = {
                               'g':trans_db.joinpath('g_HSC.res'),
                               'r':trans_db.joinpath('r_HSC.res'),
                               'i':trans_db.joinpath('i_HSC.res'),
                               'z':trans_db.joinpath('z_HSC.res'),
                               'Y':trans_db.joinpath('Y_HSC.res'),
                              }
                              
    instrument_name = 'HSC'
    
    color_dict = {
                  'g':'C2',
                  'r':'C3',
                  'i':'C7',
                  'z':'C5',
                  'Y':'C6',
                 }
    
    pass

class SDSSTransmission(TransmissionBase):
    
    transmission_dictionary = {
                               'u':trans_db.joinpath('u_SDSS.res'),
                               'g':trans_db.joinpath('g_SDSS.res'),
                               'r':trans_db.joinpath('r_SDSS.res'),
                               'i':trans_db.joinpath('i_SDSS.res'),
                               'z':trans_db.joinpath('z_SDSS.res'),
                              }
                              
    instrument_name = 'SDSS'
    
    color_dict = {
                  'u':'C0',
                  'g':'C2',
                  'r':'C3',
                  'i':'C7',
                  'z':'C5',
                 }
    
    pass

# dictionary of implemented transmission systems
implemented_transmissions = {'decam':DECamTransmission, 'sdss':SDSSTransmission, 'hsc':HSCTransmission}




