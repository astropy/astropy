# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains functions for reading and writing HDF5 tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`table_io` for more details.
"""

from __future__ import print_function

import os

import numpy as np

from ... import log

HDF5_SIGNATURE = b'\x89HDF\r\n\x1a\n'

__all__ = ['read_table_hdf5', 'write_table_hdf5']


def is_hdf5(origin, path, fileobj, *args, **kwargs):

    if fileobj is not None:
        pos = fileobj.tell()
        signature = fileobj.read(8)
        fileobj.seek(pos)
        return signature == HDF5_SIGNATURE
    elif path is not None:
        return path.endswith('.hdf5') or path.endswith('.h5')

    try:
        import h5py
    except ImportError:
        return False
    else:
        return isinstance(args[0], (h5py.highlevel.File, h5py.highlevel.Group))


def read_table_hdf5(input, path=None):
    """
    Read a Table object from an HDF5 file

    This requires `h5py <http://alfven.org/wp/hdf5-for-python/>`_ to be
    installed.

    Parameters
    ----------
    input : str or `h5py.highlevel.File` or `h5py.highlevel.Group`
        If a string, the filename to read the table from. If an h5py object,
        either the file or the group object to read the table from.
    path : str
        The path from which to read the table inside the HDF5 file.
        This should be relative to the input file or group.
    """

    try:
        import h5py
    except ImportError:
        raise Exception("h5py is required to read and write HDF5 files")

    if path is None:
        raise ValueError("table path should be set via the path= argument")
    elif path.endswith('/'):
        raise ValueError("table path should end with table name, not /")

    if '/' in path:
        group, name = path.rsplit('/', 1)
    else:
        group, name = None, path

    if isinstance(input, (h5py.highlevel.File, h5py.highlevel.Group)):
        f, g = None, input
        if group:
            try:
                g = g[group]
            except KeyError:
                raise IOError("Group {0} does not exist".format(group))
    else:
        f = h5py.File(input, 'r')
        try:
            g = f[group] if group else f
        except KeyError:
            f.close()
            raise IOError("Group {0} does not exist".format(group))

    # Check whether table exists
    if name not in g:
        if f is not None:
            f.close()
        raise IOError("Table {0} does not exist".format(path))

    # Read the table from the file
    dset = g[name]

    # Create a Table object
    from ...table import Table
    table = Table(np.array(dset))

    # Read the meta-data from the file
    table.meta.update(dset.attrs)

    if f is not None:
        f.close()

    return table


def write_table_hdf5(table, output, path=None, compression=False,
               append=False, overwrite=False):
    """
    Write a Table object to an HDF5 file

    This requires `h5py <http://alfven.org/wp/hdf5-for-python/>`_ to be
    installed.

    Parameters
    ----------
    output : str or `h5py.highlevel.File` or `h5py.highlevel.Group`
        If a string, the filename to write the table to. If an h5py object,
        either the file or the group object to write the table to.
    compression : bool
        Whether to compress the table inside the HDF5 file.
    path : str
        The path to which to write the table inside the HDF5 file.
        This should be relative to the input file or group.
    append : bool
        Whether to append the table to an existing HDF5 file.
    overwrite : bool
        Whether to overwrite any existing file without warning.
    """

    try:
        import h5py
    except ImportError:
        raise Exception("h5py is required to read and write HDF5 files")

    if path is None:
        raise ValueError("table path should be set via the path= argument")
    elif path.endswith('/'):
        raise ValueError("table path should end with table name, not /")

    if '/' in path:
        group, name = path.rsplit('/', 1)
    else:
        group, name = None, path

    if isinstance(output, h5py.highlevel.File) or \
       isinstance(output, h5py.highlevel.Group):
        f, g = None, output
        if group:
            try:
                g = g[group]
            except KeyError:
                g = g.create_group(group)
    else:
        if os.path.exists(output) and not append:
            if overwrite:
                os.remove(output)
            else:
                raise IOError("File exists: {0}".format(output))

        # Open the file for appending or writing
        f = h5py.File(output, 'a' if append else 'w')

        if group:
            if append:
                if group in f.keys():
                    g = f[group]
                else:
                    g = f.create_group(group)
            else:
                g = f.create_group(group)
        else:
            g = f

    # Check whether table already exists
    if name in g:
        if f is not None:
            f.close()
        raise IOError("Table {0} already exists".format(path))

    # Write the table to the file
    dset = g.create_dataset(name, data=table._data, compression=compression)

    # Write the meta-data to the file
    for key in table.meta:
        val = table.meta[key]
        try:
            dset.attrs[key] = val
        except TypeError:
            log.warn("Attribute `{0}` of type {1} cannot be written to "
                     "HDF5 files - skipping".format(key, type(val)))

    if f is not None:
        f.close()
