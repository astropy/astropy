# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains classes and functions to standardize access to
configuration files for Astropy and affiliated packages.

.. note::
    The configuration system makes use of the 'configobj' package, which stores
    configuration in a text format like that used in the standard library
    `ConfigParser`. More information and documentation for configobj can be
    found at http://www.voidspace.org.uk/python/configobj.html.
"""

from __future__ import division
import textwrap
from ..extern.configobj import configobj, validate
from contextlib import contextmanager

__all__ = ['ConfigurationItem', 'InvalidConfigurationItemWarning',
           'ConfigurationMissingWarning', 'get_config', 'save_config',
           'reload_config']


class InvalidConfigurationItemWarning(Warning):
    """ A Warning that is issued when the configuration value specified in the
    astropy configuration file does not match the type expected for that
    configuration value.
    """


class ConfigurationMissingWarning(Warning):
    """ A Warning that is issued when the configuration directory cannot be
    accessed (usually due to a permissions problem). If this warning appears,
    configuration items will be set to their defaults rather than read from the
    configuration file, and no configuration will persist across sessions.
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
    defaultvalue
        The default value for this item. If this is a list of strings, this
        item will be interpreted as an 'options' value - this item must be one
        of those values, and the first in the list will be taken as the default
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
        element (e.g. 'astropy' if this is 'astropy.config.configuration')
        will be used to determine the name of the configuration file, while
        the remaining items determine the section. If None, the package will be
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
    '[configuration]' section of astropy.cfg (located in the directory that
    `astropy.config.paths.get_config_dir` returns), or if the option is already
    set, it will take the value from the configuration file::

        from astropy.config import ConfigurationItem

        CFG_OPTION = ConfigurationItem('cfgoption',42,module='astropy.configuration')

    If called as ``CFG_OPTION()``, this will return the value ``42``, or some
    other integer if the ``astropy.cfg`` file specifies a different value.

    If this were a file ``astropy/configuration/__init__.py``, the `module`
    option would not be necessary, as it would automatically detect the correct
    module.

    """

    #this is used to make validation faster so a Validator object doesn't
    #have to be created every time
    _validator = validate.Validator()

    def __init__(self, name, defaultvalue='', description=None, cfgtype=None,
                      module=None):
        from warnings import warn
        from ..utils import find_current_module

        if module is None:
            module = find_current_module(2)
            if module is None:
                msg1 = 'Cannot automatically determine get_config module, '
                msg2 = 'because it is not called from inside a valid module'
                raise RuntimeError(msg1 + msg2)
            else:
                module = module.__name__

        self.name = name
        self.module = module
        self.description = description

        #now determine cfgtype if it is not given
        if cfgtype is None:
            if isinstance(defaultvalue, list):
                #it is an options list
                dvstr = [str(v) for v in defaultvalue]
                cfgtype = 'option(' + ', '.join(dvstr) + ')'
                defaultvalue = dvstr[0]
            elif isinstance(defaultvalue, bool):
                cfgtype = 'boolean'
            elif isinstance(defaultvalue, int):
                cfgtype = 'integer'
            elif isinstance(defaultvalue, float):
                cfgtype = 'float'
            else:
                cfgtype = 'string'
                defaultvalue = str(defaultvalue)

        self.cfgtype = cfgtype

        self._validate_val(defaultvalue)
        self.defaultvalue = defaultvalue

        #note that the actual value is stored in the ConfigObj file for this
        #package

        #this checks the current value to make sure it's valid for the type
        #as well as updating the ConfigObj with the default value, if it's not
        #actually in the ConfigObj
        try:
            self()
        except TypeError as e:
            # make sure it's a TypeError from __call__
            if 'Configuration value not valid:' in e.args[0]:
                warn(InvalidConfigurationItemWarning(*e.args))
            else:
                raise

    def set(self, value):
        """ Sets the current value of this `ConfigurationItem`.

        This also updates the comments that give the description and type
        information.

        .. note::
            This does *not* save the value of this `ConfigurationItem` to the
            configuration file.  To do that, use `ConfigurationItem.save` or
            `save_config`.

        Parameters
        ----------
        value
            The value this item should be set to.

        Raises
        ------
        TypeError
            If the provided `value` is not valid for this `ConfigurationItem`.
        """
        try:
            value = self._validate_val(value)
        except validate.ValidateError as e:
            msg = 'Provided value for configuration item {0} not valid: {1}'
            raise TypeError(msg.format(self.name, e.args[0]))

        sec = get_config(self.module)

        sec[self.name] = value
        sec.comments[self.name] = self._generate_comments()

    @contextmanager
    def set_temp(self, value):
        """
        Sets this item to a specified value only inside a while loop.

        Use as::
            ITEM = ConfigurationItem('ITEM', 'default', 'description')

            with ITEM.set_temp('newval'):
                ... do something that wants ITEM's value to be 'newval' ...

             # ITEM is now 'default' after the with block

        Parameters
        ----------
        value
            The value to set this item to inside the with block.

        """
        initval = self()
        self.set(value)
        yield
        self.set(initval)

    def save(self, value=None):
        """ Writes a value for this `ConfigurationItem` to the relevant
        configuration file.

        This also writes updated versions of the comments that give the
        description and type information.

        .. note::
            This only saves the value of this *particular* `ConfigurationItem`.
            To save all configuration settings for this package at once, see
            `save_config`.

        Parameters
        ----------
        value
            Save this value to the configuration file. If None, the current
            value of this `ConfigurationItem` will be saved.

        Raises
        ------
        TypeError
            If the provided `value` is not valid for this `ConfigurationItem`.
        """
        try:
            value = self() if value is None else self._validate_val(value)
        except validate.ValidateError as e:
            msg = 'Provided value for configuration item {0} not valid: {1}'
            raise TypeError(msg.format(self.name, e.args[0]))

        #Now find the  ConfigObj that this is based on
        baseobj = get_config(self.module)
        secname = baseobj.name
        cobj = baseobj
        #a ConfigObj's parent is itself, so we look for the parent with that
        while cobj.parent is not cobj:
            cobj = cobj.parent

        #use the current on disk version, which will be modified with the
        #given value and type/description
        newobj = configobj.ConfigObj(cobj.filename, interpolation=False)
        if secname is not None:
            if secname not in newobj:
                newobj[secname] = {}
            newsec = newobj[secname]

        newsec[self.name] = value
        newsec.comments[self.name] = self._generate_comments()
        newobj.write()

    def reload(self):
        """ Reloads the value of this `ConfigurationItem` from the relevant
        configuration file.

        Returns
        -------
        val
            The new value loaded from the configuration file.
        """
        baseobj = get_config(self.module)
        secname = baseobj.name

        cobj = baseobj
        #a ConfigObj's parent is itself, so we look for the parent with that
        while cobj.parent is not cobj:
            cobj = cobj.parent

        newobj = configobj.ConfigObj(cobj.filename, interpolation=False)
        if secname is not None:
            newobj = newobj[secname]

        baseobj[self.name] = newobj[self.name]

    def __call__(self):
        """ Returns the value of this `ConfigurationItem`

        Returns
        -------
        val
            This item's value, with a type determined by the `cfgtype`
            attribute.

        Raises
        ------
        TypeError
            If the configuration value as stored is not this item's type.
        """

        #get the value from the relevant `configobj.ConfigObj` object
        sec = get_config(self.module)
        if self.name not in sec:
            self.set(self.defaultvalue)
        val = sec[self.name]

        try:
            return self._validate_val(val)
        except validate.ValidateError as e:
            raise TypeError('Configuration value not valid:' + e.args[0])

    def _validate_val(self, val):
        """ Validates the provided value based on cfgtype and returns the
        type-cast value

        throws the underlying configobj exception if it fails
        """
        #note that this will normally use the *class* attribute `_validator`,
        #but if some arcane reason is needed for making a special one for an
        #instance or sub-class, it will be used
        return self._validator.check(self.cfgtype, val)

    def _generate_comments(self):
        comments = []
        comments.append('')  # adds a blank line before every entry
        if self.description is not None:
            for line in textwrap.wrap(self.description, width=76):
                comments.append(line)
        if self.cfgtype.startswith('option'):
            comments.append("Options: " + self.cfgtype[7:-1])
        return comments


# this dictionary stores the master copy of the ConfigObj's for each
# root package
_cfgobjs = {}


def get_config(packageormod=None, reload=False):
    """ Gets the configuration object or section associated with a particular
    package or module.

    Parameters
    -----------
    packageormod : str or None
        The package for which to retrieve the configuration object. If a
        string, it must be a valid package name, or if None, the package from
        which this function is called will be used.

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
    from warnings import warn

    from .paths import get_config_dir
    from ..utils import find_current_module

    if packageormod is None:
        packageormod = find_current_module(2)
        if packageormod is None:
            msg1 = 'Cannot automatically determine get_config module, '
            msg2 = 'because it is not called from inside a valid module'
            raise RuntimeError(msg1 + msg2)
        else:
            packageormod = packageormod.__name__

    packageormodspl = packageormod.split('.')
    rootname = packageormodspl[0]
    secname = '.'.join(packageormodspl[1:])

    cobj = _cfgobjs.get(rootname, None)

    if cobj is None:
        try:
            cfgfn = join(get_config_dir(), rootname + '.cfg')
            _cfgobjs[rootname] = cobj = configobj.ConfigObj(cfgfn,
                interpolation=False)
        except (IOError, OSError) as e:
            msg1 = 'Configuration defaults will be used, and configuration '
            msg2 = 'cannot be saved due to '
            errstr = '' if len(e.args) < 1 else (':' + str(e.args[0]))
            wmsg = msg1 + msg2 + e.__class__.__name__ + errstr
            warn(ConfigurationMissingWarning(wmsg))

            #This caches the object, so if the file becomes acessible, this
            #function won't see it unless the module is reloaded
            _cfgobjs[rootname] = cobj = configobj.ConfigObj(interpolation=False)

    if secname:  # not the root package
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
    the 'astropy.config.configuration' module).

    .. note::
        To save only a single item, use the `ConfigurationItem.save` method -
        this will save all options in the current session that may have been
        changed.

    Parameters
    ----------
    packageormod : str or None
        The package or module name - see `get_config` for details.
    """
    sec = get_config(packageormod)
    #look for the section that is its own parent - that's the base object
    while sec.parent is not sec:
        sec = sec.parent
    sec.write()


def reload_config(packageormod=None):
    """ Reloads configuration settings from a configuration file for the root
    package of the requested package/module.

    This overwrites any changes that may have been made in `ConfigurationItem`
    objects.  This applies for any items that are based on this file, which is
    determined by the *root* package of `packageormod` (e.g. 'astropy.cfg' for
    the 'astropy.config.configuration' module).

    Parameters
    ----------
    packageormod : str or None
        The package or module name - see `get_config` for details.
    """
    sec = get_config(packageormod)
    #look for the section that is its own parent - that's the base object
    while sec.parent is not sec:
        sec = sec.parent
    sec.reload()


def get_config_items(packageormod=None):
    """ Returns the `ConfigurationItem` objects associated with a particular
    module.

    Parameters
    ----------
    packageormod : str or None
        The package or module name or None to get the current module's items.

    Returns
    -------
    configitems : dict
        A dictionary where the keys are the name of the items as the are named
        in the module, and the values are the associated `ConfigurationItem`
        objects.

    """
    import sys
    from inspect import ismodule

    from ..utils import find_current_module

    if packageormod is None:
        packageormod = find_current_module(2)
        if packageormod is None:
            msg1 = 'Cannot automatically determine get_config module, '
            msg2 = 'because it is not called from inside a valid module'
            raise RuntimeError(msg1 + msg2)
    elif isinstance(packageormod, basestring):
        __import__(packageormod)
        packageormod = sys.modules[packageormod]
    elif ismodule(packageormod):
        pass
    else:
        raise TypeError('packageormod in get_config_items is invalid')

    configitems = {}
    for n, obj in packageormod.__dict__.iteritems():
        objcls = obj.__class__
        fqn = objcls.__module__ + '.' + objcls.__name__
        if fqn == 'astropy.config.configuration.ConfigurationItem':
            configitems[n] = obj

    return configitems


def _fix_section_blank_lines(sec, recurse=True, gotoroot=True):
    """
    Adds a blank line to the comments of any sections in the requested sections,
    recursing into subsections if `recurse` is True. If `gotoroot` is True,
    this first goes to the root of the requested section, just like
    `save_config` and `reload_config` - this does nothing if `sec` is a
    configobj already.
    """

    if not hasattr(sec, 'sections'):
        sec = get_config(sec)

        #look for the section that is its own parent - that's the base object
        if gotoroot:
            while sec.parent is not sec:
                sec = sec.parent

    for isec, snm in enumerate(sec.sections):
        comm = sec.comments[snm]
        if len(comm) == 0 or comm[-1] != '':
            if sec.parent is sec and isec == 0:
                pass  # don't do it for first section
            else:
                comm.append('')
        if recurse:
            _fix_section_blank_lines(sec[snm], True, False)


def _generate_all_config_items(pkgornm=None, reset_to_default=False):
    """ Given a root package name or package, this function simply walks
    through all the subpackages and modules, which should populate any
    ConfigurationItem objects defined at the module level. If
    `reset_to_default` is True, it also sets all of the items to their default
    values, regardless of what the file's value currently is. It then saves the
    `ConfigObj`.

    If `pkgname` is None, it determines the package based on the root package
    of the function where this function is called. Be a bit cautious about
    this, though - this might not always be what you want.
    """
    from pkgutil import get_loader, walk_packages
    from types import ModuleType

    from ..utils import find_current_module

    if pkgornm is None:
        pkgornm = find_current_module(1).__name__.split('.')[0]

    if isinstance(pkgornm, basestring):
        package = get_loader(pkgornm).load_module(pkgornm)
    elif isinstance(pkgornm, ModuleType) and '__init__' in pkgornm.__file__:
        package = pkgornm
    else:
        msg = '_generate_all_config_items was not given a package/package name'
        raise TypeError(msg)

    for imper, nm, ispkg in walk_packages(package.__path__,
                                          package.__name__ + '.'):
        mod = imper.find_module(nm)
        if reset_to_default:
            for v in mod.__dict__.itervalues():
                if isinstance(v, ConfigurationItem):
                    v.set(v.defaultvalue)

    _fix_section_blank_lines(package.__name__, True, True)
    save_config(package.__name__)
