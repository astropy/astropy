# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
import os
import re

from .base import IORegistryError

__all__ = ['UnifiedReadWriteMethod', 'UnifiedReadWrite']


# -----------------------------------------------------------------------------

class UnifiedReadWrite:
    """Base class for the worker object used in unified read() or write() methods.

    This lightweight object is created for each `read()` or `write()` call
    via ``read`` / ``write`` descriptors on the data object class.  The key
    driver is to allow complete format-specific documentation of available
    method options via a ``help()`` method, e.g. ``Table.read.help('fits')``.

    Subclasses must define a ``__call__`` method which is what actually gets
    called when the data object ``read()`` or ``write()`` method is called.

    For the canonical example see the `~astropy.table.Table` class
    implementation (in particular the ``connect.py`` module there).

    Parameters
    ----------
    instance : object
        Descriptor calling instance or None if no instance
    cls : type
        Descriptor calling class (either owner class or instance class)
    method_name : str
        Method name, e.g. 'read' or 'write'
    registry : ``_UnifiedIORegistryBase`` or None, optional
        The IO registry.
    """
    def __init__(self, instance, cls, method_name, registry=None):
        if registry is None:
            from astropy.io.registry.compat import default_registry as registry

        self._registry = registry
        self._instance = instance
        self._cls = cls
        self._method_name = method_name  # 'read' or 'write'

    @property
    def registry(self):
        """Unified I/O registry instance."""
        return self._registry

    def help(self, format=None, out=None):
        """Output help documentation for the specified unified I/O ``format``.

        By default the help output is printed to the console via ``pydoc.pager``.
        Instead one can supplied a file handle object as ``out`` and the output
        will be written to that handle.

        Parameters
        ----------
        format : str
            Unified I/O format name, e.g. 'fits' or 'ascii.ecsv'
        out : None or path-like
            Output destination (default is stdout via a pager)
        """
        cls = self._cls
        method_name = self._method_name

        # Get reader or writer function associated with the registry
        get_func = (self._registry.get_reader if method_name == 'read'
                    else self._registry.get_writer)
        try:
            if format:
                read_write_func = get_func(format, cls)
        except IORegistryError as err:
            reader_doc = 'ERROR: ' + str(err)
        else:
            if format:
                # Format-specific
                header = ("{}.{}(format='{}') documentation\n"
                          .format(cls.__name__, method_name, format))
                doc = read_write_func.__doc__
            else:
                # General docs
                header = f'{cls.__name__}.{method_name} general documentation\n'
                doc = getattr(cls, method_name).__doc__

            reader_doc = re.sub('.', '=', header)
            reader_doc += header
            reader_doc += re.sub('.', '=', header)
            reader_doc += os.linesep
            if doc is not None:
                reader_doc += inspect.cleandoc(doc)

        if out is None:
            import pydoc
            pydoc.pager(reader_doc)
        else:
            out.write(reader_doc)

    def list_formats(self, out=None):
        """Print a list of available formats to console (or ``out`` filehandle)

        out : None or file handle object
            Output destination (default is stdout via a pager)
        """
        tbl = self._registry.get_formats(self._cls, self._method_name.capitalize())
        del tbl['Data class']

        if out is None:
            tbl.pprint(max_lines=-1, max_width=-1)
        else:
            out.write('\n'.join(tbl.pformat(max_lines=-1, max_width=-1)))

        return out


# -----------------------------------------------------------------------------

class UnifiedReadWriteMethod(property):
    """Descriptor class for creating read() and write() methods in unified I/O.

    The canonical example is in the ``Table`` class, where the ``connect.py``
    module creates subclasses of the ``UnifiedReadWrite`` class.  These have
    custom ``__call__`` methods that do the setup work related to calling the
    registry read() or write() functions.  With this, the ``Table`` class
    defines read and write methods as follows::

      read = UnifiedReadWriteMethod(TableRead)
      write = UnifiedReadWriteMethod(TableWrite)

    Parameters
    ----------
    func : `~astropy.io.registry.UnifiedReadWrite` subclass
        Class that defines read or write functionality

    """
    # We subclass property to ensure that __set__ is defined and that,
    # therefore, we are a data descriptor, which cannot be overridden.
    # This also means we automatically inherit the __doc__ of fget (which will
    # be a UnifiedReadWrite subclass), and that this docstring gets recognized
    # and properly typeset by sphinx (which was previously an issue; see
    # gh-11554).
    # We override __get__ to pass both instance and class to UnifiedReadWrite.
    def __get__(self, instance, owner_cls):
        return self.fget(instance, owner_cls)
