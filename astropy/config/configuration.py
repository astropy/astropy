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

import hashlib
import io
import textwrap

from contextlib import contextmanager
from os import path
from warnings import warn

from ..extern.configobj import configobj, validate
from ..utils.exceptions import AstropyWarning
from .paths import get_config_dir


__all__ = ['ConfigurationItem', 'InvalidConfigurationItemWarning',
           'ConfigurationMissingWarning', 'get_config', 'save_config',
           'reload_config']


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
    `astropy.config.get_config_dir` returns), or if the option is already
    set, it will take the value from the configuration file::

        from astropy.config import ConfigurationItem

        CFG_OPTION = ConfigurationItem('cfgoption',42,module='astropy.configuration')

    If called as ``CFG_OPTION()``, this will return the value ``42``, or some
    other integer if the ``astropy.cfg`` file specifies a different value.

    If this were a file ``astropy/configuration/__init__.py``, the `module`
    option would not be necessary, as it would automatically detect the correct
    module.

    """

    # this is used to make validation faster so a Validator object doesn't
    # have to be created every time
    _validator = validate.Validator()

    def __init__(self, name, defaultvalue='', description=None, cfgtype=None,
                 module=None):
        from ..utils import find_current_module
        from ..utils import isiterable

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
            else:
                cfgtype = 'string'
                defaultvalue = str(defaultvalue)

        self.cfgtype = cfgtype

        self._validate_val(defaultvalue)
        self.defaultvalue = defaultvalue

        # note that the actual value is stored in the ConfigObj file for this
        # package

        # this checks the current value to make sure it's valid for the type
        # as well as updating the ConfigObj with the default value, if it's not
        # actually in the ConfigObj
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
        try:
            yield
        finally:
            self.set(initval)

    def save(self, value=None):
        """
        Removed in astropy 0.4.
        """
        raise NotImplementedError(
            "The ability to save config options was removed in astropy 0.4. "
            "To change config settings, edit '{0}' directly.".
            format(get_config_filename(self.module)))

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
        # a ConfigObj's parent is itself, so we look for the parent with that
        while cobj.parent is not cobj:
            cobj = cobj.parent

        newobj = configobj.ConfigObj(cobj.filename, interpolation=False)
        if secname is not None:
            newobj = newobj[secname]

        baseobj[self.name] = newobj[self.name]

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

        # get the value from the relevant `configobj.ConfigObj` object
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
        # note that this will normally use the *class* attribute `_validator`,
        # but if some arcane reason is needed for making a special one for an
        # instance or sub-class, it will be used
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


def get_config_filename(packageormod=None):
    """
    Get the filename of the config file associated with the given
    package or module.
    """
    cfg = get_config(packageormod)
    while cfg.parent is not cfg:
        cfg = cfg.parent
    return cfg.filename


def get_config(packageormod=None):
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
        if _ASTROPY_SETUP_:
            # There's no reason to use anything but the default config
            cobj = configobj.ConfigObj(interpolation=False)
        else:
            try:
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
