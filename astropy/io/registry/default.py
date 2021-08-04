# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools

from .core import UnifiedIORegistry

__all__ = ['register_reader', 'register_writer', 'register_identifier',
           'identify_format', 'get_reader', 'get_writer', 'read', 'write',
           'get_formats', 'delay_doc_updates',
           # added to scope
           'unregister_reader', 'unregister_writer', 'unregister_identifier']

__doctest_skip__ = ['register_identifier']

default_registry = UnifiedIORegistry()

_identifiers = default_registry._identifiers
_readers = default_registry._readers
_writers = default_registry._writers


def _make_io_func(func_name):

    @functools.wraps(getattr(default_registry, func_name))
    def wrapper(*args, registry=None, **kwargs):
        if registry is None:
            registry = default_registry
        method = getattr(registry, func_name)
        return method(*args, **kwargs)

    return wrapper


register_identifier = _make_io_func("register_identifier")
unregister_identifier = _make_io_func("unregister_identifier")
identify_format = _make_io_func("identify_format")
register_reader = _make_io_func("register_reader")
unregister_reader = _make_io_func("unregister_reader")
get_reader = _make_io_func("get_reader")
read = _make_io_func("read")
register_writer = _make_io_func("register_writer")
unregister_writer = _make_io_func("unregister_writer")
get_writer = _make_io_func("get_writer")
write = _make_io_func("write")
get_formats = _make_io_func("get_formats")
delay_doc_updates = _make_io_func("delay_doc_updates")
