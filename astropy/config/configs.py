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
from ..extern.configobj import configobj,validate

__all__ = ['ConfigurationItem','InvalidConfigurationItemWarning','get_config',
           'save_config','reload_config']


class InvalidConfigurationItemWarning(Warning):
    """ A Warning that is issued when the configuration value specified in the
    astropy configuration file does not match the type expected for that 
    configuration value.
    """

class ConfigurationItem(object):
    """ A setting and associated value stored in the astropy configuration 
    files.
    
    These objects are typically defined at the top of astropy subpackages
    or affiliated packages, and store values or option settings that can be 
    modified by the user to 
    
    Parameters
    ----------
    name : str
        The (case-sensitive) name of this parameter, as shown in the 
        configuration file.
    defaultvalue : 
        The default value for this item. If this is a list of strings, this item
        will be interpreted as an 'options' value - this item must be one of
        those values, and the first in the list will be taken as the default
        value.
    description : str or None
        A description of this item (will be shown as a comment in the 
        configuration file)
    cfgtype : str or None
        A type specifier like those used as the *values* of a particular key in
        a `configspec` file of `configobj`. If None, the type will be inferred
        from the default value.
    module : str or None
        The full module name that this item is associated with. The first
        element (e.g. 'astropy' if this is 'astropy.config.configs') will be
        used to determine the name of the configuration file, while the
        remaining items determine the section. If None, the package will be
        inferred from the package within whiich this object's initializer is
        called.
        
    Raises
    ------
    RuntimeError
        If `module` is None, but the module this item is created from cannot 
        be determined.
        
    Examples
    --------
    The following example will create an item 'cfgoption = 42' in the 
    '[configs]' section of astropy.cfg (located in the directory that
    `astropy.config.paths.get_config_dir` returns), or if the option is already
    set, it will take the value from the configuration file::
        
        from astropy.config import ConfigurationItem
        
        CFG_OPTION = ConfigurationItem('cfgoption',42,module='astropy.configs')
        
    If called as ``CFG_OPTION()``, this will return the value ``42``, or some
    other integer if the ``astropy.cfg`` file specifies a different value.    
    
    If this were in the ``astropy/configs/__init__.py`` file, the `module` 
    option would not be necessary, as it would automatically detect the correct
    module.
    
    """
    
    #this is used to make validation faster so a Validator object doesn't
    #have to be created every time
    _validator = validate.Validator()
    
    def __init__(self,name,defaultvalue='',description=None,cfgtype=None,
                      module=None):
        
        if module is None:
            module = _find_current_module(2)
            if module is None:
                msg1 = 'Cannot automatically determine get_config module, '
                msg2 = 'because it is not called from inside a valid module'
                raise RuntimeError(msg1+msg2)
            else:
                module = module.__name__
        self.module = module
        raise NotImplementedError

    def set(self,value):
        """ Sets the current value of this `ConfigurationItem`.
        
        .. note::
            This does *not* save the value of this `ConfigurationItem` to the
            configuration file.  To do that, use `ConfigurationItem.save` or
            `save_config`.
        
        Parameter
        ---------
        value : 
            The value this item should be set to. If 
            
        Raises
        ------
        TypeError
            If the provided `value` is not valid for this `ConfigurationItem`.
        """
        raise NotImplementedError
    
    def save(self,value=None):
        """ Writes a value for this `ConfigurationItem` to the relevant
        configuration file.
        
        .. note::
            This only saves the value of this *particular* `ConfigurationItem`.
            To save all configuration settings for this package at once, see
            `save_config`.
        
        Parameters
        ----------
        value : 
            Save this value to the configuration file. If None, the current
            value of this `ConfigurationItem` will be saved.
        """
        raise NotImplementedError
    
    def reload(self):
        """ Reloads the value of this `ConfigurationItem` from the relevant 
        configuration file. 
        
        Returns
        -------
        val :
            The new value loaded from the configuration file.
        """
        raise NotImplementedError
    
    def __call__(self):
        """ Returns the value of this `ConfigurationItem`
        
        Returns
        val : 
            This item's value, with a type determined by the `cfgtype` 
            attribute.
        """
        return self._validate_val(self._retrieve_val())
    
    def _retrieve_val(self):
        """ Gets and returns the value from the `configobj.ConfigObj` object
        for this item
        """
        raise NotImplementedError
        
    def _validate_val(self,val):
        """ Validates the provided value based on cfgtype and returns the 
        type-cast value
        """
        #note that this will normally use the *class* attribute `_validator`,
        #but if some arcane reason is needed for making a special one for an
        #instance or sub-class, it will be used
        return self._validator.check(self.cfgtype,val)
        
        
# this dictionary stores the master copy of the ConfigObj's for each 
# root package
_cfgobjs = {}
def get_config(packageormod=None,reload=False):
    """ Gets the configuration object or section associated with a particular
    package or module.
    
    Parameters
    -----------
    packageormod : str or None
        The package for which to retrieve the configuration object. If a string,
        it must be a valid package name, or if None, the package from which this
        function is called will be used.
        
    Returns
    -------
    cfgobj : `configobj.ConfigObj` or `configobj.Section`
        If the requested package is a base package, this will be the
        `configobj.ConfigObj` for that package, or if it is a subpackage or
        module, it will return the relevant `configobj.Section` object.
        
    Raises
    ------
    RuntimeError
        If `package` is None, but the package this item is created from cannot 
        be determined.
    
    """
    from os.path import join
    from .paths import get_config_dir
    
    if packageormod is None:
        packageormod = _find_current_module(2)
        if packageormod is None:
            msg1 = 'Cannot automatically determine get_config module, '
            msg2 = 'because it is not called from inside a valid module'
            raise RuntimeError(msg1 + msg2)
        else:
            packageormod = packageormod.__name__
    
    packageormodspl = packageormod.split('.')
    rootname = packageormodspl[0]
    secname = '.'.join(packageormodspl[1:])
    
    cobj = _cfgobjs.get(rootname,None)
    if cobj is None:
        cfgfn = join(get_config_dir(),rootname+'.cfg')
        _cfgobjs[rootname] = cobj = configobj.ConfigObj(cfgfn)
    
    if secname: #not the root package
        if secname not in cobj:
            cobj[secname] = {}
        return cobj[secname]
    else:
        return cobj

def save_config(packageormod=None):
    """ Saves all configuration settings to the configuration file for the
    root package of the requested package/module.
    
    This overwrites any configuration items that have been changed in
    `ConfigurationItem` objects that are based on the configuration file
    determined by the *root* package of `packageormod` (e.g. 'astropy.cfg' for
    the 'astropy.config.configs' module).
    
    .. note::
        To save only a single item, use the `ConfigurationItem.save` method - 
        this will save all options in the current session that may have been 
        changed.
    
    Parameters
    ----------
    packageormod : str or None
        The package or module name - see `get_config` for details. 
    """
    get_config(packageormod).write()
    
def reload_config(packageormod=None):
    """ Reloads configuration settings from a configuration file for the root
    package of the requested package/module.
    
    This overwrites any changes that may have been made in `ConfigurationItem` 
    objects.  This applies for any items that are based on this file, which is
    determined by the *root* package of `packageormod` (e.g. 'astropy.cfg' for 
    the 'astropy.config.configs' module).
    
    Parameters
    ----------
    packageormod : str or None
        The package or module name - see `get_config` for details. 
    """
    get_config(packageormod).reload()
    
def _find_current_module(depth=1):
    """ Determines the module this function is called from.  `depth` specifies
    how far back to go in the call stack - e.g. 1 indicates the package of this 
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
    
def _generate_all_config_items(package=None,reset_to_default=False):
    """ Given a root package or package name, this function simple walks through
    all the subpackages and modules, which should populate any ConfigurationItem
    objects defined at the module level. IF `reset_to_Default` is True, it also
    sets all of the items to their default values, regardless of what the file's
    value currently is. It then saves the `ConfigObj`.
    """
    import pkgutil
    
    if isinstance(package,basestring):
        package = pkgutil.find_module(package).load_module(package) 
    
    for imper,nm,ispkg in pkgutil.walk_packages(package.__path__,package.__name__+'.'):
        mod = imper.load_module(nm)
        if reset_to_default:
            for v in mod.__dict__.values():
                if isinstance(v,ConfigurationItem):
                    v.set(v.defaultvalue)
    save_config(package.__name__)
        