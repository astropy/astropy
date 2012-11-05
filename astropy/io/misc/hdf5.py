# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os

import numpy as np

from ... import log

HDF5_SIGNATURE = '\x89HDF\r\n\x1a\n'


def is_hdf5(origin, args, kwargs):

    try:
        import h5py
    except ImportError:
        raise Exception("h5py is required to read and write HDF5 files")

    if isinstance(args[0], h5py.highlevel.File) or \
        isinstance(args[0], h5py.highlevel.Group):
        return True
    elif isinstance(args[0], basestring):
        if os.path.exists(args[0]):
            with open(args[0], 'rb') as f:
                if f.read(8) == HDF5_SIGNATURE:
                    return True
                else:
                    return False
        elif args[0].endswith('.hdf5'):
            return True

    return False


def read_hdf5(input, name=None, group=""):
    """
    Read a Table object from an HDF5 file

    Parameters
    ----------
    input : str or h5py.highlevel.File or h5py.highlevel.Group
        If a string, the filename to read the table from. If an h5py object,
        either the file or the group object to read the table from.
    group : str
        The group to read the table from inside the HDF5 file. This can
        only be used if the ``input`` argument is a string.
    name : str
        The table name in the file.
    """

    try:
        import h5py
    except ImportError:
        raise Exception("h5py is required to read and write HDF5 files")

    if name is None:
        raise ValueError("table name should be set via the name= argument")
    elif '/' in name:
        raise ValueError("table name should not contain any '/'")

    if isinstance(input, h5py.highlevel.File) or \
       isinstance(input, h5py.highlevel.Group):
        f, g = None, input
        if group:
            try:
                g = g[group]
            except KeyError:
                raise Exception("Group {0} does not exist".format(group))
    else:
        f = h5py.File(input, 'r')
        g = f[group] if group else f

    # Check whether table exists
    if name not in g.keys():
        raise Exception("Table {0}/{1} does not exist".format(group, name))

    # Read the table from the file
    dset = g[name]

    # Create a Table object
    from ...table import Table
    table = Table(np.array(dset))

    # Read the meta-data from the file
    for key in dset.attrs:
        table.meta[key] = dset.attrs[key]

    if f is not None:
        f.close()

    return table


def write_hdf5(table, output, name=None, compression=False, group="",
               append=False, overwrite=False):
    """
    Write a Table object to an HDF5 file

    Parameters
    ----------
    output : str or h5py.highlevel.File or h5py.highlevel.Group
        If a string, the filename to write the table to. If an h5py object,
        either the file or the group object to write the table to.
    name : str
        The table name in the file.
    compression : bool
        Whether to compress the table inside the HDF5 file.
    group : str
        The group to write the table to inside the HDF5 file, relative
        to the output (i.e. if a h5py group object is passed in
        `output`, then `group` is the path relative to the group.)
    append : bool
        Whether to append the table to an existing HDF5 file.
    overwrite : bool
        Whether to overwrite any existing file without warning.
    """

    try:
        import h5py
    except ImportError:
        raise Exception("h5py is required to read and write HDF5 files")

    if name is None:
        raise ValueError("table name should be set via the name= argument")
    elif '/' in name:
        raise ValueError("table name should not contain any '/'")

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
                raise Exception("File exists: {0}".format(output))

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
    if name in g.keys():
        raise Exception("Table {0}/{1} already exists".format(group, name))

    # Write the table to the file
    dset = g.create_dataset(name, data=table._data, compression=compression)

    # Write the meta-data to the file
    for key in table.meta:

        if isinstance(table.meta[key], basestring):
            # Use np.string_ to ensure that fixed-length attributes are used.
            dset.attrs[key] = np.string_(table.meta[key])
        else:
            try:
                dset.attrs[key] = table.meta[key]
            except TypeError:
                log.warn("Attribute `{0}` of type {1} cannot be written to HDF5 files - skipping".format(key, type(table.meta[key])))

    if f is not None:
        f.close()

from astropy.table.io_registry import register_reader, register_writer, register_identifier

register_reader('hdf5', read_hdf5)
register_writer('hdf5', write_hdf5)
register_identifier('hdf5', is_hdf5)
