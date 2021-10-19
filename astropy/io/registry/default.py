# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .core import UnifiedIORegistry

__all__ = ['register_reader', 'register_writer', 'register_identifier',
           'identify_format', 'get_reader', 'get_writer', 'read', 'write',
           'get_formats', 'delay_doc_updates',
           # added to scope
           'unregister_reader', 'unregister_writer', 'unregister_identifier']

__doctest_skip__ = ['register_identifier']

default_registry = UnifiedIORegistry()

# identify
_identifiers = default_registry._identifiers
register_identifier = default_registry.register_identifier
unregister_identifier = default_registry.unregister_identifier
identify_format = default_registry.identify_format

# read
_readers = default_registry._readers
register_reader = default_registry.register_reader
unregister_reader = default_registry.unregister_reader
get_reader = default_registry.get_reader
read = default_registry.read

# write
_writers = default_registry._writers
register_writer = default_registry.register_writer
unregister_writer = default_registry.unregister_writer
get_writer = default_registry.get_writer
write = default_registry.write

# utils
get_formats = default_registry.get_formats
delay_doc_updates = default_registry.delay_doc_updates
