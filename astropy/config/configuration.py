# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains classes and functions to standardize access to
configuration files for Astropy and affiliated packages.

.. note::
    The configuration system makes use of the 'configobj' package, which stores
    configuration in a text format like that used in the standard library
    `ConfigParser`. More information and documentation for configobj can be
    found at http://www.voidspace.org.uk/python/configobj.html.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six

from contextlib import contextmanager
import hashlib
import io
from os import path
import sys
from warnings import warn

from ..extern.configobj import configobj, validate
from ..utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from ..utils import find_current_module
from ..utils.misc import InheritDocstrings
from .paths import get_config_dir


__all__ = ['ConfigurationItem', 'InvalidConfigurationItemWarning',
           'ConfigurationMissingWarning', 'get_config', 'save_config',
           'reload_config', 'ConfigNamespace', 'ConfigItem', 'ConfigAlias']


class InvalidConfigurationItemWarning(AstropyWarning):
    """ A Warning that is issued when the configuration value specified in the
    astropy configuration file does not match the type expected for that
    configuration value.
    """


class ConfigurationMissingWarning(AstropyWarning):
    """ A Warning that is issued when the configuration directory cannot be
    accessed (usually due to a permissions problem). If this warning appears,
    configuration items will be set to their defaults rather than read from the
    configuration file, and no configuration will persist across sessions.
    """


# these are not in __all__ because it's not intended that a user ever see them
class ConfigurationDefaultMissingError(ValueError):
    """ An exception that is raised when the configuration defaults (which
    should be generated at build-time) are missing.
    """


# this is used in astropy/__init__.py
class ConfigurationDefaultMissingWarning(AstropyWarning):
    """ A warning that is issued when the configuration defaults (which
    should be generated at build-time) are missing.
    """


class ConfigurationChangedWarning(AstropyWarning):
    """
    A warning that the configuration options have changed.
    """


class _ConfigNamespaceMeta(type):
    def __init__(cls, name, bases, dict):
        if cls.__bases__[0] is object:
            return

        for key, val in six.iteritems(dict):
            if isinstance(val, ConfigItem):
                val.name = key


@six.add_metaclass(_ConfigNamespaceMeta)
class ConfigNamespace(object):
    """
    A namespace of configuration items.  Each subpackage with
    configuration items should define a subclass of this class,
    containing `ConfigItem` instances as members.

    For example::

        class Conf(_config.ConfigNamespace):
            unicode_output = _config.ConfigItem(
                False,
                'Use Unicode characters when outputting values, ...')
            use_color = _config.ConfigItem(
                sys.platform != 'win32',
                'When True, use ANSI color escape sequences when ...',
                aliases=['astropy.utils.console.USE_COLOR'])
        conf = Conf()
    """
    def set_temp(self, attr, value):
        """
        Temporarily set a configuration value.

        Parameters
        ----------
        attr : str
            Configuration item name

        value : object
            The value to set temporarily.

        Examples
        --------
        >>> import astropy
        >>> with astropy.conf.set_temp('use_color', False):
        ...     pass
        ...     # console output will not contain color
        >>> # console output contains color again...
        """
        if hasattr(self, attr):
            return self.__class__.__dict__[attr].set_temp(value)
        raise AttributeError("No configuration parameter '{0}'".format(attr))

    def reload(self, attr=None):
        """
        Reload a configuration item from the configuration file.

        Parameters
        ----------
        attr : str, optional
            The name of the configuration parameter to reload.  If not
            provided, reload all configuration parameters.
        """
        if attr is not None:
            if hasattr(self, attr):
                return self.__class__.__dict__[attr].reload()
            raise AttributeError("No configuration parameter '{0}'".format(attr))

        for item in six.itervalues(self.__class__.__dict__):
            if isinstance(item, ConfigItem):
                item.reload()

    def reset(self, attr=None):
        """
        Reset a configuration item to its default.

        Parameters
        ----------
        attr : str, optional
            The name of the configuration parameter to reload.  If not
            provided, reset all configuration parameters.
        """
        if attr is not None:
            if hasattr(self, attr):
                prop = self.__class__.__dict__[attr]
                prop.set(prop.defaultvalue)
                return
            raise AttributeError("No configuration parameter '{0}'".format(attr))

        for item in six.itervalues(self.__class__.__dict__):
            if isinstance(item, ConfigItem):
                item.set(item.defaultvalue)


@six.add_metaclass(InheritDocstrings)
class ConfigItem(object):
    """
    A setting and associated value stored in a configuration file.

    These objects should be created as members of
    `ConfigNamespace` subclasses, for example::

        class _Conf(config.ConfigNamespace):
            unicode_output = config.ConfigItem(
                'unicode_output', False,
                'Use Unicode characters when outputting values, and writing widgets '
                'to the console.')
        conf = _Conf()

    Parameters
    ----------
    defaultvalue : object, optional
        The default value for this item. If this is a list of strings, this
        item will be interpreted as an 'options' value - this item must be one
        of those values, and the first in the list will be taken as the default
        value.

    description : str or None, optional
        A description of this item (will be shown as a comment in the
        configuration file)

    cfgtype : str or None, optional
        A type specifier like those used as the *values* of a particular key in
        a `configspec` file of `configobj`. If None, the type will be inferred
        from the default value.

    module : str or None, optional
        The full module name that this item is associated with. The first
        element (e.g. 'astropy' if this is 'astropy.config.configuration')
        will be used to determine the name of the configuration file, while
        the remaining items determine the section. If None, the package will be
        inferred from the package within whiich this object's initializer is
        called.

    aliases : str, or list of str, optional
        The deprecated location(s) of this configuration item.  If the
        config item is not found at the new location, it will be
        searched for at all of the old locations.

    Raises
    ------
    RuntimeError
        If `module` is None, but the module this item is created from cannot
        be determined.
    """

    # this is used to make validation faster so a Validator object doesn't
    # have to be created every time
    _validator = validate.Validator()

    def __init__(self, defaultvalue='', description=None, cfgtype=None,
                 module=None, aliases=None):
        from ..utils import isiterable

        if module is None:
            module = find_current_module(2)
            if module is None:
                msg1 = 'Cannot automatically determine get_config module, '
                msg2 = 'because it is not called from inside a valid module'
                raise RuntimeError(msg1 + msg2)
            else:
                module = module.__name__

        self.module = module
        self.description = description
        self.__doc__ = description

        # now determine cfgtype if it is not given
        if cfgtype is None:
            if (isiterable(defaultvalue) and not
                    isinstance(defaultvalue, six.string_types)):
                # it is an options list
                dvstr = [str(v) for v in defaultvalue]
                cfgtype = 'option(' + ', '.join(dvstr) + ')'
                defaultvalue = dvstr[0]
            elif isinstance(defaultvalue, bool):
                cfgtype = 'boolean'
            elif isinstance(defaultvalue, int):
                cfgtype = 'integer'
            elif isinstance(defaultvalue, float):
                cfgtype = 'float'
            elif isinstance(defaultvalue, six.string_types):
                cfgtype = 'string'
                defaultvalue = str(defaultvalue)

        self.cfgtype = cfgtype

        self._validate_val(defaultvalue)
        self.defaultvalue = defaultvalue

        if aliases is None:
            self.aliases = []
        elif isinstance(aliases, six.string_types):
            self.aliases = [aliases]
        else:
            self.aliases = aliases

    def __set__(self, obj, value):
        return self.set(value)

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return self()

    def set(self, value):
        """ Sets the current value of this `ConfigItem`.

        This also updates the comments that give the description and type
        information.

        Parameters
        ----------
        value
            The value this item should be set to.

        Raises
        ------
        TypeError
            If the provided `value` is not valid for this `ConfigItem`.
        """
        try:
            value = self._validate_val(value)
        except validate.ValidateError as e:
            msg = 'Provided value for configuration item {0} not valid: {1}'
            raise TypeError(msg.format(self.name, e.args[0]))

        sec = get_config(self.module)

        sec[self.name] = value

    @contextmanager
    def set_temp(self, value):
        """
        Sets this item to a specified value only inside a with block.

        Use as::
            ITEM = ConfigItem('ITEM', 'default', 'description')

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
        try:
            yield
        finally:
            self.set(initval)

    def reload(self):
        """ Reloads the value of this `ConfigItem` from the relevant
        configuration file.

        Returns
        -------
        val
            The new value loaded from the configuration file.
        """
        self.set(self.defaultvalue)
        baseobj = get_config(self.module, True)
        secname = baseobj.name

        cobj = baseobj
        # a ConfigObj's parent is itself, so we look for the parent with that
        while cobj.parent is not cobj:
            cobj = cobj.parent

        newobj = configobj.ConfigObj(cobj.filename, interpolation=False)
        if secname is not None:
            if secname not in newobj:
                return baseobj.get(self.name)
            newobj = newobj[secname]

        if self.name in newobj:
            baseobj[self.name] = newobj[self.name]
        return baseobj.get(self.name)

    def __repr__(self):
        out = '<{0}: name={1!r} value={2!r} at 0x{3:x}>'.format(
            self.__class__.__name__, self.name, self(), id(self))
        return out

    def __str__(self):
        out = '\n'.join(('{0}: {1}',
                         '  cfgtype={2!r}',
                         '  defaultvalue={3!r}',
                         '  description={4!r}',
                         '  module={5}',
                         '  value={6!r}'))
        out = out.format(self.__class__.__name__, self.name, self.cfgtype,
                         self.defaultvalue, self.description, self.module,
                         self())
        return out

    def __call__(self):
        """ Returns the value of this `ConfigItem`

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
        options = []
        sec = get_config(self.module)
        if self.name in sec:
            options.append((sec[self.name], self.module, self.name))

        for alias in self.aliases:
            module, name = alias.rsplit('.', 1)
            sec = get_config(module)
            if name in sec:
                warn(
                    "Config parameter '{0}' in section [{1}] is deprecated. "
                    "Use '{2}' in section [{3}] instead.".format(
                        name, module, self.name, self.module),
                    AstropyDeprecationWarning)
                options.append((sec[name], module, name))

        if len(options) == 0:
            self.set(self.defaultvalue)
            options.append((self.defaultvalue, None, None))

        if len(options) > 1:
            warn(
                "Config parameter '{0}' in section [{1}] is defined by "
                "more than one alias: {2}".format(
                    self.name, self.module,
                    ['.'.join(x[1], x[2]) for x in options]))

        val = options[0][0]

        try:
            return self._validate_val(val)
        except validate.ValidateError as e:
            raise TypeError('Configuration value not valid:' + e.args[0])

    def _validate_val(self, val):
        """ Validates the provided value based on cfgtype and returns the
        type-cast value

        throws the underlying configobj exception if it fails
        """
        # note that this will normally use the *class* attribute `_validator`,
        # but if some arcane reason is needed for making a special one for an
        # instance or sub-class, it will be used
        return self._validator.check(self.cfgtype, val)


class ConfigurationItem(ConfigItem):
    """
    A backward-compatibility layer to support the old
    `ConfigurationItem` API.  The only difference between this and
    `ConfigItem` is that this requires an explicit name to be set as
    the first argument.
    """
    # REMOVE in astropy 0.5

    def __init__(self, name, defaultvalue='', description=None, cfgtype=None,
                 module=None, aliases=None):
        warn(
            "ConfigurationItem has been deprecated in astropy 0.4. "
            "Use ConfigItem objects as members of ConfigNamespace subclasses "
            "instead.  See ConfigNamespace for an example.",
            AstropyDeprecationWarning)

        # We have to do the automatic module determination here, not
        # just in ConfigItem, otherwise the extra stack frame will
        # make it come up with the wrong answer.
        if module is None:
            module = find_current_module(2)
            if module is None:
                msg1 = 'Cannot automatically determine get_config module, '
                msg2 = 'because it is not called from inside a valid module'
                raise RuntimeError(msg1 + msg2)
            else:
                module = module.__name__

        super(ConfigurationItem, self).__init__(
            defaultvalue=defaultvalue,
            description=description,
            cfgtype=cfgtype,
            module=module,
            aliases=aliases)
        self.name = name

    def save(self, value=None):
        """
        Removed in astropy 0.4.
        """
        raise NotImplementedError(
            "The ability to save config options was removed in astropy 0.4. "
            "To change config settings, edit '{0}' directly.".
            format(get_config_filename(self.module)))


class ConfigAlias(ConfigItem):
    """
    A class that exists to support backward compatibility only.

    This is an alias for a `ConfigItem` that has been moved elsewhere.
    It inherits from `ConfigItem` only because it implements the same
    interface, not because any of the methods are reused.

    Parameters
    ----------
    since : str
        The version in which the configuration item was moved.

    old_name : str
        The old name of the configuration item.  This should be the
        name of the variable in Python, not in the configuration file.

    new_name : str
        The new name of the configuration item.  This is both the name
        of the item in Python and in the configuration file (since as of
        astropy 0.4, those are always the same thing).

    old_module : str, optional
        A fully-qualified, dot-separated path to the module in which
        the configuration item used to be defined.  If not provided, it
        is the name of the module in which `ConfigAlias` is called.

    new_module : str, optional
        A fully-qualified, dot-separated path to the module in which
        the configuration item is now defined.  If not provided, it is
        the name of the module in which `ConfigAlias` is called.  This
        string should not contain the `.conf` object.  For example, if
        the new configuration item is in `astropy.conf.use_unicode`, this
        value only needs to be `astropy`.
    """
    # REMOVE in astropy 0.5

    def __init__(self, since, old_name, new_name, old_module=None, new_module=None):
        if old_module is None:
            old_module = find_current_module(2)
            if old_module is None:
                msg1 = 'Cannot automatically determine get_config module, '
                msg2 = 'because it is not called from inside a valid module'
                raise RuntimeError(msg1 + msg2)
            else:
                old_module = old_module.__name__

        if new_module is None:
            new_module = old_module

        self._since = since
        self._old_name = old_name
        self._new_name = new_name
        self._old_module = old_module
        self._new_module = new_module

    def _deprecation_warning(self):
        warn(
            "Since {0}, config parameter '{1}.{2}' is deprecated. "
            "Use '{3}.conf.{4}' instead.".format(
                self._since,
                self._old_module, self._old_name,
                self._new_module, self._new_name),
            AstropyDeprecationWarning)

    def _get_target(self):
        if self._new_module not in sys.modules:
            __import__(self._new_module)
        mod = sys.modules[self._new_module]
        cfg = getattr(mod, 'conf')
        return cfg

    def set(self, value):
        self._deprecation_warning()
        setattr(self._get_target(), self._new_name, value)

    def set_temp(self, value):
        self._deprecation_warning()
        return self._get_target().set_temp(self._new_name, value)

    def save(self, value=None):
        self._deprecation_warning()
        return self._get_target().save(value)

    def reload(self):
        self._deprecation_warning()
        return self._get_target().reload(self._new_name)

    def __repr__(self):
        return repr(getattr(self._get_target().__class__, self._new_name))

    def __str__(self):
        return str(getattr(self._get_target().__class__, self._new_name))

    def __call__(self):
        self._deprecation_warning()
        return getattr(self._get_target(), self._new_name)


# this dictionary stores the master copy of the ConfigObj's for each
# root package
_cfgobjs = {}


def get_config_filename(packageormod=None):
    """
    Get the filename of the config file associated with the given
    package or module.
    """
    cfg = get_config(packageormod)
    while cfg.parent is not cfg:
        cfg = cfg.parent
    return cfg.filename


# This is used by testing to override the config file, so we can test
# with various config files that exercise different features of the
# config system.
_override_config_file = None


def get_config(packageormod=None, reload=False):
    """ Gets the configuration object or section associated with a particular
    package or module.

    Parameters
    -----------
    packageormod : str or None
        The package for which to retrieve the configuration object. If a
        string, it must be a valid package name, or if None, the package from
        which this function is called will be used.

    reload : bool, optional
        Reload the file, even if we have it cached.

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

    if cobj is None or reload:
        if _ASTROPY_SETUP_:
            # There's no reason to use anything but the default config
            cobj = configobj.ConfigObj(interpolation=False)
        else:
            try:
                # This feature is intended only for use by the unit tests
                if _override_config_file is not None:
                    cfgfn = _override_config_file
                else:
                    cfgfn = path.join(get_config_dir(), rootname + '.cfg')
                cobj = configobj.ConfigObj(cfgfn, interpolation=False)
            except (IOError, OSError) as e:
                msg = ('Configuration defaults will be used due to ')
                errstr = '' if len(e.args) < 1 else (':' + str(e.args[0]))
                msg += e.__class__.__name__ + errstr
                warn(ConfigurationMissingWarning(msg))

                # This caches the object, so if the file becomes accessible, this
                # function won't see it unless the module is reloaded
                cobj = configobj.ConfigObj(interpolation=False)

        _cfgobjs[rootname] = cobj

    if secname:  # not the root package
        if secname not in cobj:
            cobj[secname] = {}
        return cobj[secname]
    else:
        return cobj


def save_config(packageormod=None, filename=None):
    """
    Removed in astropy 0.4.
    """
    raise NotImplementedError(
        "The ability to save config options was removed in astropy 0.4. "
        "To change config settings, edit '{0}' directly.".
        format(get_config_filename(packageormod)))


def reload_config(packageormod=None):
    """ Reloads configuration settings from a configuration file for the root
    package of the requested package/module.

    This overwrites any changes that may have been made in `ConfigItem`
    objects.  This applies for any items that are based on this file, which is
    determined by the *root* package of `packageormod` (e.g. 'astropy.cfg' for
    the 'astropy.config.configuration' module).

    Parameters
    ----------
    packageormod : str or None
        The package or module name - see `get_config` for details.
    """
    sec = get_config(packageormod, True)
    # look for the section that is its own parent - that's the base object
    while sec.parent is not sec:
        sec = sec.parent
    sec.reload()


def is_unedited_config_file(filename):
    """
    Determines if a config file can be safely replaced because it doesn't
    actually contain any meaningful content.

    To meet this criteria, the config file must be either:

    - All comments or completely empty

    - An exact match to a "legacy" version of the config file prior to
      Astropy 0.4, when APE3 was implemented and the config file
      contained commented-out values by default.
    """
    with open(filename, 'rb') as fd:
        content = fd.read()

    # First determine if the config file has any effective content
    buffer = io.BytesIO(content)
    buffer.seek(0)
    raw_cfg = configobj.ConfigObj(buffer, interpolation=True)
    for v in six.itervalues(raw_cfg):
        if len(v):
            break
    else:
        return True

    # Now determine if it matches the md5sum of a known, unedited
    # config file.
    # TODO: How does this work with Windows line endings...?  Probably
    # doesn't...
    known_configs = set([
        '5df7e409425e5bfe7ed041513fda3288',  # v0.3
        '8355f99a01b3bdfd8761ef45d5d8b7e5',  # v0.2
        '4ea5a84de146dc3fcea2a5b93735e634'   # v0.2.1, v0.2.2, v0.2.3, v0.2.4, v0.2.5
    ])

    md5 = hashlib.md5()
    md5.update(content)
    digest = md5.hexdigest()
    return digest in known_configs


# this is not in __all__ because it's not intended that a user uses it
def update_default_config(pkg, default_cfg_dir_or_fn, version=None):
    """
    Checks if the configuration file for the specified package exists,
    and if not, copy over the default configuration.  If the
    configuration file looks like it has already been edited, we do
    not write over it, but instead write a file alongside it named
    ``pkg.version.cfg`` as a "template" for the user.

    Parameters
    ----------
    pkg : str
        The package to be updated.
    default_cfg_dir_or_fn : str
        The filename or directory name where the default configuration file is.
        If a directory name, `pkg`.cfg will be used in that directory.
    version : str, optional
        The current version of the given package.  If not provided, it will
        be obtained from ``pkg.__version__``.

    Returns
    -------
    updated : bool
        If the profile was updated, True, otherwise False.

    Raises
    ------
    ConfigurationDefaultMissingError
        If the default configuration could not be found.

    """
    cfgfn = get_config(pkg).filename

    if path.exists(cfgfn):
        doupdate = is_unedited_config_file(cfgfn)
    else:
        doupdate = True

    if version is None:
        mod = __import__(pkg)
        if not hasattr(mod, '__version__'):
            raise ConfigurationDefaultMissingError(
                'Could not determine version of package {0}'.format(pkg))
        version = mod.__version__

    # Don't install template files for dev versions, or we'll end up
    # spamming `~/.astropy/config`.
    if not 'dev' in version:
        template_path = path.join(
            get_config_dir(), '{0}.{1}.cfg'.format(pkg, version))
        needs_template = not path.exists(template_path)
    else:
        needs_template = False

    if doupdate or needs_template:
        if path.isdir(default_cfg_dir_or_fn):
            default_cfgfn = path.join(default_cfg_dir_or_fn, pkg + '.cfg')
        else:
            default_cfgfn = default_cfg_dir_or_fn

        if not path.isfile(default_cfgfn):
            # TODO: Since this file is in the repository now, it seems very
            # unlikely it would be missing...  Remove this?
            raise ConfigurationDefaultMissingError(
                'Requested default configuration file {0} is '
                'not a file.'.format(default_cfgfn))

        with open(default_cfgfn, 'r') as fr:
            content = fr.read()

        if needs_template:
            with open(template_path, 'w') as fw:
                fw.write(content)
            # If we just installed a new template file and we can't
            # update the main configuration file because it has user
            # changes, display a warning.
            if not doupdate:
                warn(
                    "The configuration options in {0} {1} may have changed, "
                    "your configuration file was not updated in order to "
                    "preserve local changes.  A new configuration template "
                    "has been saved to '{2}'.".format(
                        pkg, version, template_path),
                    ConfigurationChangedWarning)

        if doupdate:
            with open(cfgfn, 'w') as fw:
                fw.write(content)
            return True

    return False
