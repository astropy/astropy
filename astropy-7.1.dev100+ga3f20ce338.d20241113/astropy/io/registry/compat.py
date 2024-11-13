# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools

from .core import UnifiedIORegistry

__all__ = [  # noqa: F822
    "register_reader",
    "register_writer",
    "register_identifier",
    "unregister_reader",
    "unregister_writer",
    "unregister_identifier",
    "get_reader",
    "get_writer",
    "get_formats",
    "read",
    "write",
    "identify_format",
    "delay_doc_updates",
]

# make a default global-state registry  (not publicly scoped, but often accessed)
# this is for backward compatibility when ``io.registry`` was a file.
default_registry = UnifiedIORegistry()
# also need to expose the enclosed registries
_identifiers = default_registry._identifiers
_readers = default_registry._readers
_writers = default_registry._writers


def _make_io_func(method_name):
    """Makes a function for a method on UnifiedIORegistry.

    .. todo::

        Make kwarg "registry" not hidden.

    Returns
    -------
    wrapper : callable
        Signature matches method on UnifiedIORegistry.
        Accepts (hidden) kwarg "registry". default is ``default_registry``.
    """

    @functools.wraps(getattr(default_registry, method_name))
    def wrapper(*args, registry=None, **kwargs):
        # written this way in case ever controlled by ScienceState
        if registry is None:
            registry = default_registry
        # get and call bound method from registry instance
        return getattr(registry, method_name)(*args, **kwargs)

    return wrapper


# =============================================================================
# JIT function creation and lookup (PEP 562)


def __dir__():
    dir_out = list(globals())
    return sorted(dir_out + __all__)


def __getattr__(method: str):
    if method in __all__:
        return _make_io_func(method)

    raise AttributeError(f"module {__name__!r} has no attribute {method!r}")
