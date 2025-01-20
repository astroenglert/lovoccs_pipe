
import configparser

from importlib import resources as impresources
from pathlib import Path

__config_directory__ = impresources.files(__package__).joinpath('config_files')

# a config function for loading filepaths
def filepath_configs():
    '''
    
    A function for loading a dictionary of keyed filepaths from a config file
    
    Args:
        None
    
    Returns:
        config: Config; a config object storing filepaths, can be interacted w. as if a dictionary
    
    '''
    
    # load the filepath
    filepath_config = __config_directory__.joinpath('filepath_configs.ini')
    
    # read the file
    config = configparser.ConfigParser()
    config.read(filepath_config)
    
    return config
    

# a config function for loading filepaths
def pipeline_configs():
    '''
    
    A function for loading a dictionary of keyed filepaths from a config file
    
    Args:
        None
    
    Returns:
        config: Config; a config object storing options, can be interacted w. as if a dictionary
    
    '''
    
    # load the filepath
    filepath_config = __config_directory__.joinpath('pipeline_configs.ini')
    
    # read the file
    config = configparser.ConfigParser()
    config.read(filepath_config)
    
    return config


