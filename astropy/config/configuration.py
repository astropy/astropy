# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains classes and functions to standardize access to
configuration files for Astropy and affiliated packages.

.. note::
    The configuration system makes use of the 'configobj' package, which stores
    configuration in a text format like that used in the standard library
    `ConfigParser`. More information and documentation for configobj can be
    found at https://configobj.readthedocs.io .
"""

from __future__ import annotations

import contextlib
import importlib
import io
import operator
import os
import pkgutil
import warnings
from contextlib import contextmanager, nullcontext
from functools import reduce
from inspect import getdoc
from pathlib import Path
from textwrap import TextWrapper
from typing import TYPE_CHECKING
from warnings import warn

import numpy as np

from astropy.extern.configobj import configobj, validate
from astropy.utils import find_current_module, silence
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning

from .paths import get_config_dir_path

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Final

__all__ = (
    "ConfigItem",
    "ConfigNamespace",
    "InvalidConfigurationItemWarning",
    "create_config_file",
    "generate_config",
    "get_config",
    "reload_config",
)


class InvalidConfigurationItemWarning(AstropyWarning):
    """A Warning that is issued when the configuration value specified in the
    astropy configuration file does not match the type expected for that
    configuration value.
    """


# these are not in __all__ because it's not intended that a user ever see them
class ConfigurationDefaultMissingError(ValueError):
    """An exception that is raised when the configuration defaults (which
    should be generated at build-time) are missing.
    """


# this is used in astropy/__init__.py
class ConfigurationDefaultMissingWarning(AstropyWarning):
    """A warning that is issued when the configuration defaults (which
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

    def __iter__(self) -> Generator[str, None, None]:
        for key, val in self.__class__.__dict__.items():
            if isinstance(val, ConfigItem):
                yield key

    def __str__(self):
        if (docstring := getdoc(self)) is not None:
            header = f"{docstring}\n\n"
        else:
            current_module = str(find_current_module(2)).split("'")[1]
            header = f"Configuration parameters for `{current_module}`\n\n"
        return header + "\n\n".join(map(str, self.values()))

    keys = __iter__
    """Iterate over configuration item names."""

    def values(self) -> Generator[ConfigItem, None, None]:
        """Iterate over configuration item values."""
        for val in self.__class__.__dict__.values():
            if isinstance(val, ConfigItem):
                yield val

    def items(self) -> Generator[tuple[str, ConfigItem], None, None]:
        """Iterate over configuration item ``(name, value)`` pairs."""
        for key, val in self.__class__.__dict__.items():
            if isinstance(val, ConfigItem):
                yield key, val

    def help(self, name: str | None = None) -> None:
        """Print info about configuration items.

        Parameters
        ----------
        name : `str`, optional
            Name of the configuration item to be described. If no name is
            provided then info about all the configuration items will be
            printed.

        Examples
        --------
        >>> from astropy import conf
        >>> conf.help("unicode_output")
        ConfigItem: unicode_output
          cfgtype='boolean'
          defaultvalue=False
          description='When True, use Unicode characters when outputting values, and displaying widgets at the console.'
          module=astropy
          value=False
        """
        if name is None:
            print(self)
        else:
            try:
                print(type(self).__dict__[name])
            except KeyError:
                raise KeyError(
                    f"'{name}' is not among configuration items {tuple(self)}"
                ) from None

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

    def reset(self, attr: str | None = None) -> None:
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

    rootname = "astropy"
    """
    Rootname sets the base path for all config files.
    """

    def __init__(
        self, defaultvalue="", description=None, cfgtype=None, module=None, aliases=None
    ):
        if module is None:
            module = find_current_module(2)
            if module is None:
                msg1 = "Cannot automatically determine get_config module, "
                msg2 = "because it is not called from inside a valid module"
                raise RuntimeError(msg1 + msg2)
            else:
                module = module.__name__

        self.module = module
        self.description = description
        self.__doc__ = description

        # now determine cfgtype if it is not given
        if cfgtype is None:
            if np.iterable(defaultvalue) and not isinstance(defaultvalue, str):
                # it is an options list
                dvstr = [str(v) for v in defaultvalue]
                cfgtype = "option(" + ", ".join(dvstr) + ")"
                defaultvalue = dvstr[0]
            elif isinstance(defaultvalue, bool):
                cfgtype = "boolean"
            elif isinstance(defaultvalue, int):
                cfgtype = "integer"
            elif isinstance(defaultvalue, float):
                cfgtype = "float"
            elif isinstance(defaultvalue, str):
                cfgtype = "string"
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
            raise TypeError(
                f"Provided value for configuration item {self.name} not valid:"
                f" {e.args[0]}"
            )

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
        """Reloads the value of this ``ConfigItem`` from the relevant
        configuration file.

        Returns
        -------
        val : object
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

    def __repr__(self) -> str:
        return (
            f"<{self.__class__.__name__}: name={self.name!r} value={self()!r} at"
            f" 0x{id(self):x}>"
        )

    def __str__(self) -> str:
        return "\n".join(
            (
                f"{self.__class__.__name__}: {self.name}",
                f"  cfgtype={self.cfgtype!r}",
                f"  defaultvalue={self.defaultvalue!r}",
                f"  description={self.description!r}",
                f"  module={self.module}",
                f"  value={self()!r}",
            )
        )

    def __call__(self):
        """Returns the value of this ``ConfigItem``.

        Returns
        -------
        val : object
            This item's value, with a type determined by the ``cfgtype``
            attribute.

        Raises
        ------
        TypeError
            If the configuration value as stored is not this item's type.

        """

        def section_name(section):
            if section == "":
                return "at the top-level"
            else:
                return f"in section [{section}]"

        options = []
        sec = get_config(self.module, rootname=self.rootname)
        if self.name in sec:
            options.append((sec[self.name], self.module, self.name))

        for alias in self.aliases:
            module, name = alias.rsplit(".", 1)
            sec = get_config(module, rootname=self.rootname)
            if "." in module:
                filename, module = module.split(".", 1)
            else:
                filename = module
                module = ""
            if name in sec:
                if "." in self.module:
                    new_module = self.module.split(".", 1)[1]
                else:
                    new_module = ""
                warn(
                    f"Config parameter '{name}' {section_name(module)} of the file"
                    f" '{get_config_filename(filename, rootname=self.rootname)}' is"
                    f" deprecated. Use '{self.name}'"
                    f" {section_name(new_module)} instead.",
                    AstropyDeprecationWarning,
                )
                options.append((sec[name], module, name))

        if len(options) == 0:
            self.set(self.defaultvalue)
            options.append((self.defaultvalue, None, None))

        if len(options) > 1:
            filename, sec = self.module.split(".", 1)
            warn(
                f"Config parameter '{self.name}' {section_name(sec)} of the file"
                f" '{get_config_filename(filename, rootname=self.rootname)}' is given"
                " by more than one alias"
                f" ({', '.join(['.'.join(x[1:3]) for x in options if x[1] is not None])})."
                " Using the first.",
                AstropyDeprecationWarning,
            )

        val = options[0][0]

        try:
            return self._validate_val(val)
        except validate.ValidateError as e:
            raise TypeError(f"Configuration value not valid: {e.args[0]}")

    def _validate_val(self, val):
        """Validates the provided value based on cfgtype and returns the
        type-cast value.

        throws the underlying configobj exception if it fails
        """
        # note that this will normally use the *class* attribute `_validator`,
        # but if some arcane reason is needed for making a special one for an
        # instance or sub-class, it will be used
        return self._validator.check(self.cfgtype, val)


# this dictionary stores the primary copy of the ConfigObj's for each
# root package
_cfgobjs: Final[dict[str, configobj.ConfigObj]] = {}


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
    """Gets the configuration object or section associated with a particular
    package or module.

    Parameters
    ----------
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
            msg1 = "Cannot automatically determine get_config module, "
            msg2 = "because it is not called from inside a valid module"
            raise RuntimeError(msg1 + msg2)
        else:
            packageormod = packageormod.__name__

        _autopkg = True

    else:
        _autopkg = False

    packageormodspl = packageormod.split(".")
    pkgname = packageormodspl[0]
    secname = ".".join(packageormodspl[1:])

    if rootname is None:
        if _autopkg:
            rootname = pkgname
        else:
            rootname = "astropy"  # so we don't break affiliated packages

    cobj = _cfgobjs.get(pkgname)

    if cobj is None or reload:
        cfgfn = None
        try:
            # This feature is intended only for use by the unit tests
            if _override_config_file is not None:
                cfgfn = Path(_override_config_file)
            else:
                cfgfn = (
                    get_config_dir_path(rootname=rootname)
                    .joinpath(pkgname)
                    .with_suffix(".cfg")
                )
            cobj = configobj.ConfigObj(str(cfgfn), interpolation=False)
        except OSError:
            # This can happen when HOME is not set
            cobj = configobj.ConfigObj(interpolation=False)

        # This caches the object, so if the file becomes accessible, this
        # function won't see it unless the module is reloaded
        _cfgobjs[pkgname] = cobj

    if secname:  # not the root package
        if secname not in cobj:
            cobj[secname] = {}
        return cobj[secname]
    else:
        return cobj


def _recursive_subclasses(class_):
    """
    Return all subclasses of all subclasses.
    """
    subclasses = class_.__subclasses__()
    if not subclasses:
        return []
    next = reduce(operator.concat, [_recursive_subclasses(cls) for cls in subclasses])
    return [*next, *subclasses]


def generate_config(pkgname="astropy", filename=None, verbose=False):
    """Generates a configuration file, from the list of `ConfigItem`
    objects for each subpackage.

    .. versionadded:: 4.1

    Parameters
    ----------
    pkgname : str or None
        The package for which to retrieve the configuration object.
    filename : str or file-like or None
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
        warnings.simplefilter("ignore", category=filter_warnings)
        for mod in pkgutil.walk_packages(
            path=package.__path__, prefix=package.__name__ + "."
        ):
            if mod.module_finder.path.endswith(("test", "tests")) or mod.name.endswith(
                "setup_package"
            ):
                # Skip test and setup_package modules
                continue
            if mod.name.split(".")[-1].startswith("_"):
                # Skip private modules
                continue

            with contextlib.suppress(ImportError):
                importlib.import_module(mod.name)

    wrapper = TextWrapper(initial_indent="## ", subsequent_indent="## ", width=78)

    if filename is None:
        filename = get_config_filename(pkgname)

    with contextlib.ExitStack() as stack:
        if isinstance(filename, (str, os.PathLike)):
            fp = stack.enter_context(open(filename, "w"))
        else:
            # assume it's a file object, or io.StringIO
            fp = filename

        subclasses = _recursive_subclasses(ConfigNamespace)
        processed = set()

        for conf in sorted(subclasses, key=lambda x: x.__module__):
            mod = conf.__module__

            # Skip modules for other packages, e.g. astropy modules that
            # would be imported when running the function for astroquery.
            if mod.split(".")[0] != pkgname:
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
                        modname = item.module.replace(f"{pkgname}.", "")
                        fp.write(f"[{modname}]\n\n")
                    print_module = False

                fp.write(wrapper.fill(item.description) + "\n")
                if isinstance(item.defaultvalue, (tuple, list)):
                    if len(item.defaultvalue) == 0:
                        fp.write(f"# {item.name} = ,\n\n")
                    elif len(item.defaultvalue) == 1:
                        fp.write(f"# {item.name} = {item.defaultvalue[0]},\n\n")
                    else:
                        fp.write(
                            f"# {item.name} ="
                            f" {','.join(map(str, item.defaultvalue))}\n\n"
                        )
                else:
                    fp.write(f"# {item.name} = {item.defaultvalue}\n\n")


def reload_config(packageormod=None, rootname=None):
    """Reloads configuration settings from a configuration file for the root
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
    actually contain any meaningful content, i.e. if it contains only comments
    or is completely empty.
    """
    buffer = io.StringIO(content)
    raw_cfg = configobj.ConfigObj(buffer, interpolation=True)
    # If any of the items is set, return False
    return not any(len(v) > 0 for v in raw_cfg.values())


def create_config_file(pkg, rootname="astropy", overwrite=False):
    """
    Create the default configuration file for the specified package.
    If the file already exists, it is updated only if it has not been
    modified.  Otherwise the ``overwrite`` flag is needed to overwrite it.

    Parameters
    ----------
    pkg : str
        The package to be updated.
    rootname : str
        Name of the root configuration directory.
    overwrite : bool
        Force updating the file if it already exists.

    Returns
    -------
    updated : bool
        If the profile was updated, `True`, otherwise `False`.

    """
    # local import to prevent using the logger before it is configured
    from astropy.logger import log

    cfgfn = Path(get_config_filename(pkg, rootname=rootname))

    # generate the default config template
    template_content = io.StringIO()
    generate_config(pkg, template_content)
    template_content.seek(0)
    template_content = template_content.read()

    doupdate = True

    # if the file already exists, check that it has not been modified
    if cfgfn is not None and cfgfn.is_file():
        with open(cfgfn, encoding="latin-1") as fd:
            content = fd.read()

        doupdate = is_unedited_config_file(content, template_content)

    if doupdate or overwrite:
        with open(cfgfn, "w", encoding="latin-1") as fw:
            fw.write(template_content)
        log.info(f"The configuration file has been successfully written to {cfgfn}")
        return True
    elif not doupdate:
        log.warning(
            "The configuration file already exists and seems to "
            "have been customized, so it has not been updated. "
            "Use overwrite=True if you really want to update it."
        )

    return False
