# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains functions for reading and writing HDF5 tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`table_io` for more details.
"""

from __future__ import print_function

import os
import warnings

import numpy as np
from ...utils.exceptions import AstropyUserWarning

HDF5_SIGNATURE = b'\x89HDF\r\n\x1a\n'

__all__ = ['read_table_hdf5', 'write_table_hdf5']


def _find_all_structured_arrays(handle):
    """
    Find all sturctured arrays in an HDF5 file
    """
    import h5py
    structured_arrays = []

    def append_structured_arrays(name, obj):
        if isinstance(obj, h5py.Dataset) and obj.dtype.kind == 'V':
            structured_arrays.append(name)
    handle.visititems(append_structured_arrays)
    return structured_arrays


def is_hdf5(origin, filepath, fileobj, *args, **kwargs):

    if fileobj is not None:
        loc = fileobj.tell()
        try:
            signature = fileobj.read(8)
        finally:
            fileobj.seek(loc)
        return signature == HDF5_SIGNATURE
    elif filepath is not None:
        return filepath.endswith('.hdf5') or filepath.endswith('.h5')

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
    installed. If more than one table is present in the HDF5 file or group, the
    first table is read in and a warning is displayed.

    Parameters
    ----------
    input : str or `h5py.highlevel.File` or `h5py.highlevel.Group` or `h5py.highlevel.Dataset`
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

    if isinstance(input, (h5py.highlevel.File, h5py.highlevel.Group)):

        # If a path was specified, follow the path

        if path:
            try:
                input = input[path]
            except KeyError:
                raise IOError("Path {0} does not exist".format(path))

        # `input` is now either a group or a dataset. If it is a group, we
        # will search for all structured arrays inside the group, and if there
        # is one we can proceed otherwise an error is raised. If it is a
        # dataset, we just proceed with the reading.

        if isinstance(input, h5py.highlevel.Group):

            # Find all structured arrays in group
            arrays = _find_all_structured_arrays(input)

            if len(arrays) == 0:
                raise ValueError("no table found in HDF5 group {0}".format(path))
            elif len(arrays) > 0:
                path = arrays[0] if path is None else path + '/' + arrays[0]
                warnings.warn("path= was not specified but multiple tables"
                              " are present, reading in first available"
                              " table (path={0})".format(path), AstropyUserWarning)
                return read_table_hdf5(input, path=path)

    elif not isinstance(input, h5py.highlevel.Dataset):

        # If a file object was passed, then we need to extract the filename
        # because h5py cannot properly read in file objects.

        if hasattr(input, 'read'):
            try:
                input = input.name
            except AttributeError:
                raise TypeError("h5py can only open regular files")

        # Open the file for reading, and recursively call read_table_hdf5 with
        # the file object and the path.

        f = h5py.File(input, 'r')

        try:
            return read_table_hdf5(f, path=path)
        finally:
            f.close()

    # If we are here, `input` should be a Dataset object, which we can now
    # convert to a Table.

    # Create a Table object
    from ...table import Table
    table = Table(np.array(input))

    # Read the meta-data from the file
    table.meta.update(input.attrs)

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

    if isinstance(output, (h5py.highlevel.File, h5py.highlevel.Group)):

        if group:
            try:
                output_group = output[group]
            except KeyError:
                output_group = output.create_group(group)
        else:
            output_group = output

    elif isinstance(output, basestring):

        if os.path.exists(output) and not append:
            if overwrite:
                os.remove(output)
            else:
                raise IOError("File exists: {0}".format(output))

        # Open the file for appending or writing
        f = h5py.File(output, 'a' if append else 'w')

        # Recursively call the write function
        try:
            return write_table_hdf5(table, f, path=path,
                                    compression=compression, append=append,
                                    overwrite=overwrite)
        finally:
            f.close()

    else:

        raise TypeError('output should be a string or an h5py File or Group object')

    # Check whether table already exists
    if name in output_group:
        raise IOError("Table {0} already exists".format(path))

    # Write the table to the file
    dset = output_group.create_dataset(name, data=table._data, compression=compression)

    # Write the meta-data to the file
    for key in table.meta:
        val = table.meta[key]
        try:
            dset.attrs[key] = val
        except TypeError:
            warnings.warn("Attribute `{0}` of type {1} cannot be written to "
                          "HDF5 files - skipping".format(key, type(val)), AstropyUserWarning)
