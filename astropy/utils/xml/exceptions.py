# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Generic XML exception and warning helpers."""

# STDLIB
import warnings

# ASTROPY
from astropy.utils.exceptions import AstropyWarning

__all__ = [
    "XMLWarning",
    "raise_exception",
    "reraise",
    "warn",
    "warn_or_raise",
    "warn_unknown_attrs",
]


def _format_message(message, name, config=None, pos=None):
    if config is None:
        config = {}
    if pos is None:
        pos = ("?", "?")
    filename = config.get("filename", "?")
    return f"{filename}:{pos[0]}:{pos[1]}: {name}: {message}"


def _suppressed_warning(warning, config, max_warnings=10, stacklevel=2):
    warning_class = type(warning)
    config.setdefault("_warning_counts", {}).setdefault(warning_class, 0)
    config["_warning_counts"][warning_class] += 1
    message_count = config["_warning_counts"][warning_class]
    if message_count <= max_warnings:
        if message_count == max_warnings:
            warning.formatted_message += (
                " (suppressing further warnings of this type...)"
            )
        warnings.warn(warning, stacklevel=stacklevel + 1)


def warn_or_raise(
    warning_class,
    exception_class=None,
    args=(),
    config=None,
    pos=None,
    stacklevel=1,
    max_warnings=10,
):
    """Warn or raise an exception depending on the verify setting."""
    if config is None:
        config = {}
    config_value = config.get("verify", "warn")
    if config_value == "exception":
        if exception_class is None:
            exception_class = warning_class
        raise_exception(exception_class, args, config, pos)
    elif config_value == "warn":
        warnings.warn(
            warning_class,
            args,
            config,
            pos,
            stacklevel=stacklevel + 1,
            max_warnings=max_warnings,
        )


def raise_exception(exception_class, args=(), config=None, pos=None):
    """Raise an exception, with proper position information if available."""
    if config is None:
        config = {}
    raise exception_class(args, config, pos)


def reraise(exc, config=None, pos=None, additional=""):
    """Raise an exception, preserving traceback and attaching position info."""
    if config is None:
        config = {}
    message = _format_message(str(exc), exc.__class__.__name__, config, pos)
    if message.split()[0] == str(exc).split()[0]:
        message = str(exc)
    if len(additional):
        message += " " + additional
    exc.args = (message,)
    raise exc


def warn(warning_class, args=(), config=None, pos=None, stacklevel=1, max_warnings=10):
    """Warn, with proper position information if available."""
    if config is None:
        config = {}
    if config.get("verify", "warn") != "ignore":
        warning = warning_class(args, config, pos)
        _suppressed_warning(
            warning, config, max_warnings=max_warnings, stacklevel=stacklevel + 1
        )


def warn_unknown_attrs(
    warning_class, element, attrs, config, pos, good_attr=None, stacklevel=1
):
    if good_attr is None:
        good_attr = []
    for attr in attrs:
        if attr not in good_attr:
            warnings.warn(
                warning_class(args=(attr, element), config=config, pos=pos),
                stacklevel=stacklevel + 1,
            )


class XMLWarning(AstropyWarning):
    """Base class for XML warnings and exceptions."""

    default_args = ()
    message_template = ""

    def __init__(self, args, config=None, pos=None):
        if config is None:
            config = {}
        if not isinstance(args, tuple):
            args = (args,)
        msg = self.message_template.format(*args)

        self.formatted_message = _format_message(
            msg, self.__class__.__name__, config, pos
        )
        Warning.__init__(self, self.formatted_message)

    def __str__(self):
        return self.formatted_message

    @classmethod
    def get_short_name(cls):
        if len(cls.default_args):
            return cls.message_template.format(*cls.default_args)
        return cls.message_template
