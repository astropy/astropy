# Licensed under a 3-clause BSD style license - see LICENSE.rst  
"""This module contains classes and functions to standardize access to 
configuration files for Astropy and affiliated packages.

.. note::
    The configuration system makes use of the 'configobj' pacakge, which stores
    configuration in a text format like that used in the standard library
    `ConfigParser`. More information and documentation for confobj can be found
    at `http://www.voidspace.org.uk/python/configobj.html`_.
"""

from __future__ import division
from ..extern.configobj import configobj

__all__ = ['ConfigurationItem','get_config']

class ConfigurationItem(object):
    """ Fill in docstring
    
    Parameters
    ----------
    name : str
        The (case-sensitive) name of this parameter, as shown in the 
        configuration file.
    defaultvalue : 
        The default value for this item.
    description : str or None
        A description of this item (will be shown as a comment in the 
        configuration file)
    package : str or None
        The full package name that this item is associated with.  The first
        element of the pacakge name (e.g. `astropy`) will be used to determine
        the name of the configuration file, while the remaining items determine
        the section.  If None, the package will be inferred from the package 
        within whiich this object's initializer is called.
        
    Raises
    ------
    RuntimeError
        If `package` is None, but the package this item is created from cannot 
        be determined.
    
    """
    
    def __init__(self,name,defaultvalue='',description=None,package=None):
        
        if package is None:
            package = _find_current_module(2)
            if package is None:
                raise RuntimeError('get_config not called from inside a valid package')
            else:
                package = package.__name__
        self.package = package
        raise NotImplementedError

    def set(self,value):
        """ Sets the current valur of this `ConfigurationItem` 
        """
        raise NotImplementedError
    
    def write(self,value=None):
        """ Saves this `ConfigurationItem` to the relevant configuration file.
        
        Parameters
        ----------
        value : 
            Save this `ConfigurationItem` with the provided value.  If None, 
            the current value will be saved to the configuration file.  
        """
        raise NotImplementedError
    
def get_config(self,package=None):
    """
    Gets the configuration item associated with a particular package.
    
    Parameters
    -----------
    package : str or None
        The package for which to retrieve the configuration object. If a string,
        it must be a valid package name, or if None, the package from which this
        function is called will be used.
        
    Returns
    -------
    cfgobj : `configobj.ConfigObj` or `configobj.Section`
        If the requested package is the base package for 
        
    Raises
    ------
    RuntimeError
        If `package` is None, but the package this item is created from cannot 
        be determined.
    
    """
    if package is None:
        package = _find_current_module(2)
        if package is None:
            raise RuntimeError('get_config not called from inside a valid package')
    raise NotImplementedError
    
def _find_current_module(depth=1):
    """
    Determines the module this function is called from.  `depth` specifies how
    far back to go in the call stack - e.g. 1 indicates the package of this 
    function's caller, 2 retreives the caller's caller, etc.
    
    Returns the module object or None if the package cannot be found.  The name
    of the module is available as the ``__name__`` attribute of the returned
    object (if it isn't None)
    """
    import inspect
    
    frm = inspect.currentframe()
    for i in range(depth):
        frm = frm.f_back
        if frm is None:
            return None
    return inspect.getmodule(frm)
    
def _generate_all_config_items(package=None):
    """
    Given a root package, this function looks for all `ConfigurationItems`
    in it and it's subpackage/modules, and saves them out to the relevant file
    with default values.  Note that this *will* overwrite any existing values.
    """
    raise NotImplementedError