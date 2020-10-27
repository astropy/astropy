# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains classes and functions to standardize access to
configuration files for Astropy and affiliated packages.

.. note::
    The configuration system makes use of the 'configobj' package, which stores
    configuration in a text format like that used in the standard library
    `ConfigParser`. More information and documentation for configobj can be
    found at http://www.voidspace.org.uk/python/configobj.html.
"""

import io
import re
import hashlib
import pathlib
import pkgutil
import warnings
import importlib
import contextlib
from os import path
from textwrap import TextWrapper
from warnings import warn
from contextlib import contextmanager

from astropy.extern.configobj import configobj, validate
from astropy.utils import find_current_module, silence
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning
from astropy.utils.introspection import resolve_name
from astropy.utils.compat.context import nullcontext

from .paths import get_config_dir

__all__ = ['InvalidConfigurationItemWarning',
           'ConfigurationMissingWarning', 'get_config',
           'reload_config', 'ConfigNamespace', 'ConfigItem',
           'generate_config']


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

        for key, val in dict.items():
            if isinstance(val, ConfigItem):
                val.name = key


class ConfigNamespace(metaclass=_ConfigNamespaceMeta):
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
    def __iter__(self):
        for key, val in self.__class__.__dict__.items():
            if isinstance(val, ConfigItem):
                yield key

    keys = __iter__
    """Iterate over configuration item names."""

    def values(self):
        """Iterate over configuration item values."""
        for val in self.__class__.__dict__.values():
            if isinstance(val, ConfigItem):
                yield val

    def items(self):
        """Iterate over configuration item ``(name, value)`` pairs."""
        for key, val in self.__class__.__dict__.items():
            if isinstance(val, ConfigItem):
                yield key, val

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
        raise AttributeError(f"No configuration parameter '{attr}'")

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
            raise AttributeError(f"No configuration parameter '{attr}'")

        for item in self.values():
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
            raise AttributeError(f"No configuration parameter '{attr}'")

        for item in self.values():
            item.set(item.defaultvalue)


class ConfigItem:
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
        inferred from the package within which this object's initializer is
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

    rootname = 'astropy'
    """
    Rootname sets the base path for all config files.
    """

    def __init__(self, defaultvalue='', description=None, cfgtype=None,
                 module=None, aliases=None):
        from astropy.utils import isiterable

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
                    isinstance(defaultvalue, str)):
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
            elif isinstance(defaultvalue, str):
                cfgtype = 'string'
                defaultvalue = str(defaultvalue)

        self.cfgtype = cfgtype

        self._validate_val(defaultvalue)
        self.defaultvalue = defaultvalue

        if aliases is None:
            self.aliases = []
        elif isinstance(aliases, str):
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

        sec = get_config(self.module, rootname=self.rootname)

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
        baseobj = get_config(self.module, True, rootname=self.rootname)
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
        out = '<{}: name={!r} value={!r} at 0x{:x}>'.format(
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
                return f'in section [{section}]'

        options = []
        sec = get_config(self.module, rootname=self.rootname)
        if self.name in sec:
            options.append((sec[self.name], self.module, self.name))

        for alias in self.aliases:
            module, name = alias.rsplit('.', 1)
            sec = get_config(module, rootname=self.rootname)
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
                    "Config parameter '{}' {} of the file '{}' "
                    "is deprecated. Use '{}' {} instead.".format(
                        name, section_name(module), get_config_filename(filename,
                                                                        rootname=self.rootname),
                        self.name, section_name(new_module)),
                    AstropyDeprecationWarning)
                options.append((sec[name], module, name))

        if len(options) == 0:
            self.set(self.defaultvalue)
            options.append((self.defaultvalue, None, None))

        if len(options) > 1:
            filename, sec = self.module.split('.', 1)
            warn(
                "Config parameter '{}' {} of the file '{}' is "
                "given by more than one alias ({}). Using the first.".format(
                    self.name, section_name(sec), get_config_filename(filename,
                                                                      rootname=self.rootname),
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


# this dictionary stores the master copy of the ConfigObj's for each
# root package
_cfgobjs = {}


def get_config_filename(packageormod=None, rootname=None):
    """
    Get the filename of the config file associated with the given
    package or module.
    """
    cfg = get_config(packageormod, rootname=rootname)
    while cfg.parent is not cfg:
        cfg = cfg.parent
    return cfg.filename


# This is used by testing to override the config file, so we can test
# with various config files that exercise different features of the
# config system.
_override_config_file = None


def get_config(packageormod=None, reload=False, rootname=None):
    """ Gets the configuration object or section associated with a particular
    package or module.

    Parameters
    -----------
    packageormod : str or None
        The package for which to retrieve the configuration object. If a
        string, it must be a valid package name, or if ``None``, the package from
        which this function is called will be used.

    reload : bool, optional
        Reload the file, even if we have it cached.

    rootname : str or None
        Name of the root configuration directory. If ``None`` and
        ``packageormod`` is ``None``, this defaults to be the name of
        the package from which this function is called. If ``None`` and
        ``packageormod`` is not ``None``, this defaults to ``astropy``.

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

        _autopkg = True

    else:
        _autopkg = False

    packageormodspl = packageormod.split('.')
    pkgname = packageormodspl[0]
    secname = '.'.join(packageormodspl[1:])

    if rootname is None:
        if _autopkg:
            rootname = pkgname
        else:
            rootname = 'astropy'  # so we don't break affiliated packages

    cobj = _cfgobjs.get(pkgname, None)

    if cobj is None or reload:
        cfgfn = None
        try:
            # This feature is intended only for use by the unit tests
            if _override_config_file is not None:
                cfgfn = _override_config_file
            else:
                cfgfn = path.join(get_config_dir(rootname=rootname), pkgname + '.cfg')
            cobj = configobj.ConfigObj(cfgfn, interpolation=False)
        except OSError as e:
            msg = ('Configuration defaults will be used due to ')
            errstr = '' if len(e.args) < 1 else (':' + str(e.args[0]))
            msg += e.__class__.__name__ + errstr
            msg += f' on {cfgfn}'
            warn(ConfigurationMissingWarning(msg))

            # This caches the object, so if the file becomes accessible, this
            # function won't see it unless the module is reloaded
            cobj = configobj.ConfigObj(interpolation=False)

        _cfgobjs[pkgname] = cobj

    if secname:  # not the root package
        if secname not in cobj:
            cobj[secname] = {}
        return cobj[secname]
    else:
        return cobj


def generate_config(pkgname='astropy', filename=None, verbose=False):
    """Generates a configuration file, from the list of `ConfigItem`
    objects for each subpackage.

    .. versionadded:: 4.1

    Parameters
    ----------
    pkgname : str or None
        The package for which to retrieve the configuration object.
    filename : str or file object or None
        If None, the default configuration path is taken from `get_config`.

    """
    if verbose:
        verbosity = nullcontext
        filter_warnings = AstropyDeprecationWarning
    else:
        verbosity = silence
        filter_warnings = Warning

    package = importlib.import_module(pkgname)
    with verbosity(), warnings.catch_warnings():
        warnings.simplefilter('ignore', category=filter_warnings)
        for mod in pkgutil.walk_packages(path=package.__path__,
                                         prefix=package.__name__ + '.'):

            if (mod.module_finder.path.endswith(('test', 'tests')) or
                    mod.name.endswith('setup_package')):
                # Skip test and setup_package modules
                continue
            if mod.name.split('.')[-1].startswith('_'):
                # Skip private modules
                continue

            with contextlib.suppress(ImportError):
                importlib.import_module(mod.name)

    wrapper = TextWrapper(initial_indent="## ", subsequent_indent='## ',
                          width=78)

    if filename is None:
        filename = get_config(pkgname).filename

    with contextlib.ExitStack() as stack:
        if isinstance(filename, (str, pathlib.Path)):
            fp = stack.enter_context(open(filename, 'w'))
        else:
            # assume it's a file object, or io.StringIO
            fp = filename

        # Parse the subclasses, ordered by their module name
        subclasses = ConfigNamespace.__subclasses__()
        processed = set()

        for conf in sorted(subclasses, key=lambda x: x.__module__):
            mod = conf.__module__

            # Skip modules for other packages, e.g. astropy modules that
            # would be imported when running the function for astroquery.
            if mod.split('.')[0] != pkgname:
                continue

            # Check that modules are not processed twice, which can happen
            # when they are imported in another module.
            if mod in processed:
                continue
            else:
                processed.add(mod)

            print_module = True
            for item in conf().values():
                if print_module:
                    # If this is the first item of the module, we print the
                    # module name, but not if this is the root package...
                    if item.module != pkgname:
                        modname = item.module.replace(f'{pkgname}.', '')
                        fp.write(f"[{modname}]\n\n")
                    print_module = False

                fp.write(wrapper.fill(item.description) + '\n')
                fp.write(f'# {item.name} = {item.defaultvalue}\n\n')


def reload_config(packageormod=None, rootname=None):
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
    rootname : str or None
        Name of the root configuration directory - see `get_config`
        for details.
    """
    sec = get_config(packageormod, True, rootname=rootname)
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
    content = re.sub(br'\njquery_url\s*=\s*[^\n]+', b'', content)

    # First determine if the config file has any effective content
    buffer = io.BytesIO(content)
    buffer.seek(0)
    raw_cfg = configobj.ConfigObj(buffer, interpolation=True)
    for v in raw_cfg.values():
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
def update_default_config(pkg, default_cfg_dir_or_fn, version=None, rootname='astropy'):
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
    rootname : str
        Name of the root configuration directory.

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

    cfgfn = get_config(pkg, rootname=rootname).filename

    with open(default_cfgfn, 'rt', encoding='latin-1') as fr:
        template_content = fr.read()

    doupdate = False
    if cfgfn is not None:
        if path.exists(cfgfn):
            with open(cfgfn, 'rt', encoding='latin-1') as fd:
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
    if version and 'dev' not in version and cfgfn is not None:
        template_path = path.join(
            get_config_dir(rootname=rootname), f'{pkg}.{version}.cfg')
        needs_template = not path.exists(template_path)
    else:
        needs_template = False

    if doupdate or needs_template:
        if needs_template:
            with open(template_path, 'wt', encoding='latin-1') as fw:
                fw.write(template_content)
            # If we just installed a new template file and we can't
            # update the main configuration file because it has user
            # changes, display a warning.
            if not identical and not doupdate:
                warn(
                    "The configuration options in {} {} may have changed, "
                    "your configuration file was not updated in order to "
                    "preserve local changes.  A new configuration template "
                    "has been saved to '{}'.".format(
                        pkg, version, template_path),
                    ConfigurationChangedWarning)

        if doupdate and not identical:
            with open(cfgfn, 'wt', encoding='latin-1') as fw:
                fw.write(template_content)
            return True

    return False
