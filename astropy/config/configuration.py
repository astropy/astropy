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
import inspect
import io
from os import path
import pkgutil
import re
import sys
import types
from warnings import warn

from ..extern.configobj import configobj, validate
from ..utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from ..utils import find_current_module
from ..utils.introspection import resolve_name
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
                False,
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
        A type specifier like those used as the *values* of a particular key
        in a ``configspec`` file of ``configobj``. If None, the type will be
        inferred from the default value.

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
        If ``module`` is `None`, but the module this item is created from
        cannot be determined.
    """

    # this is used to make validation faster so a Validator object doesn't
    # have to be created every time
    _validator = validate.Validator()
    cfgtype = None
    """
    A type specifier like those used as the *values* of a particular key in a
    ``configspec`` file of ``configobj``.
    """

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
                dvstr = [six.text_type(v) for v in defaultvalue]
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
                defaultvalue = six.text_type(defaultvalue)

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
        """
        Sets the current value of this ``ConfigItem``.

        This also updates the comments that give the description and type
        information.

        Parameters
        ----------
        value
            The value this item should be set to.

        Raises
        ------
        TypeError
            If the provided ``value`` is not valid for this ``ConfigItem``.
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
                #... do something that wants ITEM's value to be 'newval' ...
                print(ITEM)

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
        """ Reloads the value of this ``ConfigItem`` from the relevant
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
        """ Returns the value of this ``ConfigItem``

        Returns
        -------
        val
            This item's value, with a type determined by the ``cfgtype``
            attribute.

        Raises
        ------
        TypeError
            If the configuration value as stored is not this item's type.
        """
        def section_name(section):
            if section == '':
                return 'at the top-level'
            else:
                return 'in section [{0}]'.format(section)

        options = []
        sec = get_config(self.module)
        if self.name in sec:
            options.append((sec[self.name], self.module, self.name))

        for alias in self.aliases:
            module, name = alias.rsplit('.', 1)
            sec = get_config(module)
            if '.' in module:
                filename, module = module.split('.', 1)
            else:
                filename = module
                module = ''
            if name in sec:
                if '.' in self.module:
                    new_module = self.module.split('.', 1)[1]
                else:
                    new_module = ''
                warn(
                    "Config parameter '{0}' {1} of the file '{2}' "
                    "is deprecated. Use '{3}' {4} instead.".format(
                        name, section_name(module), get_config_filename(filename),
                        self.name, section_name(new_module)),
                    AstropyDeprecationWarning)
                options.append((sec[name], module, name))

        if len(options) == 0:
            self.set(self.defaultvalue)
            options.append((self.defaultvalue, None, None))

        if len(options) > 1:
            filename, sec = self.module.split('.', 1)
            warn(
                "Config parameter '{0}' {1} of the file '{2}' is "
                "given by more than one alias ({3}). Using the first.".format(
                    self.name, section_name(sec), get_config_filename(filename),
                    ', '.join([
                        '.'.join(x[1:3]) for x in options if x[1] is not None])),
                AstropyDeprecationWarning)

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
    ``ConfigItem`` is that this requires an explicit name to be set as
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
        string should not contain the ``.conf`` object.  For example, if
        the new configuration item is in ``astropy.conf.use_unicode``, this
        value only needs to be ``'astropy'``.
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
        return resolve_name(self._new_module, 'conf')

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
        string, it must be a valid package name, or if `None`, the package from
        which this function is called will be used.

    reload : bool, optional
        Reload the file, even if we have it cached.

    Returns
    -------
    cfgobj : ``configobj.ConfigObj`` or ``configobj.Section``
        If the requested package is a base package, this will be the
        ``configobj.ConfigObj`` for that package, or if it is a subpackage or
        module, it will return the relevant ``configobj.Section`` object.

    Raises
    ------
    RuntimeError
        If ``packageormod`` is `None`, but the package this item is created
        from cannot be determined.
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
            cfgfn = None
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
                msg += ' on {0}'.format(cfgfn)
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
    objects.  This applies for any items that are based on this file, which
    is determined by the *root* package of ``packageormod``
    (e.g. ``'astropy.cfg'`` for the ``'astropy.config.configuration'``
    module).

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


def is_unedited_config_file(content, template_content=None):
    """
    Determines if a config file can be safely replaced because it doesn't
    actually contain any meaningful content.

    To meet this criteria, the config file must be either:

    - All comments or completely empty

    - An exact match to a "legacy" version of the config file prior to
      Astropy 0.4, when APE3 was implemented and the config file
      contained commented-out values by default.
    """
    # We want to calculate the md5sum using universal line endings, so
    # that even if the files had their line endings converted to \r\n
    # on Windows, this will still work.

    content = content.encode('latin-1')

    # The jquery_url setting, present in 0.3.2 and later only, is
    # effectively auto-generated by the build system, so we need to
    # ignore it in the md5sum calculation for 0.3.2.
    content = re.sub(b'\njquery_url\s*=\s*[^\n]+', b'', content)

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
    known_configs = set([
        '7d4b4f1120304b286d71f205975b1286',  # v0.3.2
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
        If a directory name, ``'pkg.cfg'`` will be used in that directory.
    version : str, optional
        The current version of the given package.  If not provided, it will
        be obtained from ``pkg.__version__``.

    Returns
    -------
    updated : bool
        If the profile was updated, `True`, otherwise `False`.

    Raises
    ------
    AttributeError
        If the version number of the package could not determined.

    """

    if path.isdir(default_cfg_dir_or_fn):
        default_cfgfn = path.join(default_cfg_dir_or_fn, pkg + '.cfg')
    else:
        default_cfgfn = default_cfg_dir_or_fn

    if not path.isfile(default_cfgfn):
        # There is no template configuration file, which basically
        # means the affiliated package is not using the configuration
        # system, so just return.
        return False

    cfgfn = get_config(pkg).filename

    with io.open(default_cfgfn, 'rt', encoding='latin-1') as fr:
        template_content = fr.read()

    doupdate = False
    if cfgfn is not None:
        if path.exists(cfgfn):
            with io.open(cfgfn, 'rt', encoding='latin-1') as fd:
                content = fd.read()

            identical = (content == template_content)

            if not identical:
                doupdate = is_unedited_config_file(
                    content, template_content)
        elif path.exists(path.dirname(cfgfn)):
            doupdate = True
            identical = False

    if version is None:
        version = resolve_name(pkg, '__version__')

    # Don't install template files for dev versions, or we'll end up
    # spamming `~/.astropy/config`.
    if not 'dev' in version and cfgfn is not None:
        template_path = path.join(
            get_config_dir(), '{0}.{1}.cfg'.format(pkg, version))
        needs_template = not path.exists(template_path)
    else:
        needs_template = False

    if doupdate or needs_template:
        if needs_template:
            with io.open(template_path, 'wt', encoding='latin-1') as fw:
                fw.write(template_content)
            # If we just installed a new template file and we can't
            # update the main configuration file because it has user
            # changes, display a warning.
            if not identical and not doupdate:
                warn(
                    "The configuration options in {0} {1} may have changed, "
                    "your configuration file was not updated in order to "
                    "preserve local changes.  A new configuration template "
                    "has been saved to '{2}'.".format(
                        pkg, version, template_path),
                    ConfigurationChangedWarning)

        if doupdate and not identical:
            with io.open(cfgfn, 'wt', encoding='latin-1') as fw:
                fw.write(template_content)
            return True

    return False


# DEPRECATED FUNCTIONALITY ----------------------------------------
# Everything below this point should be removed in astropy 0.5

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

    from ..utils import find_current_module

    if packageormod is None:
        packageormod = find_current_module(2)
        if packageormod is None:
            msg1 = 'Cannot automatically determine get_config module, '
            msg2 = 'because it is not called from inside a valid module'
            raise RuntimeError(msg1 + msg2)
    elif isinstance(packageormod, six.string_types):
        __import__(packageormod)
        packageormod = sys.modules[packageormod]
    elif inspect.ismodule(packageormod):
        pass
    else:
        raise TypeError('packageormod in get_config_items is invalid')

    configitems = {}
    for n, obj in six.iteritems(packageormod.__dict__):
        # if it's not a new-style object, it's certainly not a ConfigurationItem
        if hasattr(obj, '__class__'):
            fqn = obj.__class__.__module__ + '.' + obj.__class__.__name__
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

        # look for the section that is its own parent - that's the base object
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


def _save_config(packageormod=None, filename=None):
    """ Saves all configuration settings to the configuration file for the
    root package of the requested package/module.

    This overwrites any configuration items that have been changed in
    `ConfigurationItem` objects that are based on the configuration file
    determined by the *root* package of ``packageormod`` (e.g. 'astropy.cfg'
    for the 'astropy.config.configuration' module).

    .. note::
        To save only a single item, use the `ConfigurationItem.save` method -
        this will save all options in the current session that may have been
        changed.

    Parameters
    ----------
    packageormod : str or None
        The package or module name - see `get_config` for details.

    filename : str, optional
        Save the config to a given filename instead of to the default location.

    """

    sec = get_config(packageormod)
    # look for the section that is its own parent - that's the base object
    while sec.parent is not sec:
        sec = sec.parent
    if filename is not None:
        with io.open(filename, 'w', encoding='utf-8') as f:
            sec.write(outfile=f)
    else:
        sec.write()


def generate_all_config_items(pkgornm=None, reset_to_default=False,
                              filename=None):
    """ Given a root package name or package, this function walks
    through all the subpackages and modules, which should populate any
    ConfigurationItem objects defined at the module level. If
    `reset_to_default` is True, it also sets all of the items to their default
    values, regardless of what the file's value currently is. It then saves the
    `ConfigObj`.

    Parameters
    ----------
    pkgname : str, module, or None
        The package for which to generate configuration items.  If None,
        the package of the function that calls this one will be used.

    reset_to_default : bool
        If True, the configuration items will all be set to their defaults.

    filename : str, optional
        Save the generated config items to the given filename instead of to
        the default config file path.

    Returns
    -------
    cfgfn : str
        The filename of the generated configuration item.

    """

    from ..utils import find_current_module

    unsafe_import_regex = [r'.*.setup_package']
    unsafe_import_regex = [('(' + pat + ')') for pat in _unsafe_import_regex]
    unsafe_import_regex = re.compile('|'.join(_unsafe_import_regex))

    if pkgornm is None:
        pkgornm = find_current_module(1).__name__.split('.')[0]

    if isinstance(pkgornm, six.string_types):
        package = pkgutil.get_loader(pkgornm).load_module(pkgornm)
    elif (isinstance(pkgornm, types.ModuleType) and
            '__init__' in pkgornm.__file__):
        package = pkgornm
    else:
        msg = 'generate_all_config_items was not given a package/package name'
        raise TypeError(msg)

    if hasattr(package, '__path__'):
        pkgpath = package.__path__
    elif hasattr(package, '__file__'):
        pkgpath = path.split(package.__file__)[0]
    else:
        raise AttributeError('package to generate config items for does not '
                             'have __file__ or __path__')

    prefix = package.__name__ + '.'
    for imper, nm, ispkg in pkgutil.walk_packages(pkgpath, prefix):
        if nm == 'astropy.config.tests.test_configs':
            continue
        if not unsafe_import_regex.match(nm):
            imper.find_module(nm)
            if reset_to_default:
                for cfgitem in six.itervalues(get_config_items(nm)):
                    cfgitem.set(cfgitem.defaultvalue)

    _fix_section_blank_lines(package.__name__, True, True)

    _save_config(package.__name__, filename=filename)

    if filename is None:
        return get_config(package.__name__).filename
    else:
        return filename
