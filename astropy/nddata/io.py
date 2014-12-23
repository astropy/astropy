# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the I/O mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import warnings

from ..io import registry as io_registry
from ..io import fits

from ..extern import six
from ..utils.compat.odict import OrderedDict
from ..utils.exceptions import AstropyUserWarning

__all__ = ['NDIOMixin']


class NDIOMixin(object):
    """
    Mixin class to connect NDData to the astropy input/output registry.

    This mixin adds two methods to its subclasses, ``read`` and ``write``.
    """

    def __init__(self, *args, **kwargs):
        # Take care of the regular initialization for this class.
        super(NDIOMixin, self).__init__(*args, **kwargs)

        # Register reader/writer for whatever class this is mixed in to.
        NDIOMixin._register_with_registry(self.__class__)

    @classmethod
    def _register_with_registry(cls, as_class):
        """
        Method to actually register this class with the astropy I/O registry.

        Parameters
        ----------

        as_class : class
            The class for which this I/O should be registered.
        """
        # Get Table of registered formats for this class.
        formats = io_registry.get_formats(as_class)

        # Only register if we haven't already.
        if 'fits' not in formats['Format']:
            reader = reader_for_class(as_class)
            io_registry.register_reader('fits', as_class,
                                        reader)
            io_registry.register_writer('fits', as_class,
                                        write_nddata_fits)
            io_registry.register_identifier('fits', as_class,
                                            fits.connect.is_fits)

    @classmethod
    def read(cls, *args, **kwargs):
        """
        Read and parse gridded N-dimensional data and return as an
        NDData-derived object.

        This function provides the NDDataBase interface to the astropy unified
        I/O layer.  This allows easily reading a file in the supported data
        formats.
        """
        # If the object is being constructed by a read so that _init__ has
        # never run, we may need to register read/writer.
        cls._register_with_registry(cls)

        # The underlying read function, read_nddata_fits, takes the class of
        # the object it is constructing as *its* first argument, and
        # io_registry.read removes *its* first argument before calling the
        # underlying reader, so we repeat it.
        #
        # Ignore the comment above -- that approach doesn't work because
        # the code in the i/o registry that identifies formats assumes that
        # if the first argument has a read method it must be a file-like
        # object...but it isn't. And we can't eliminate the read method,
        # obviously.
        return io_registry.read(cls, *args, **kwargs)

    def write(self, *args, **kwargs):
        """
        Write a gridded N-dimensional data object out in specified format.

        This function provides the NDDataBase interface to the astropy unified
        I/O layer.  This allows easily writing a file in the supported data
        formats.
        """
        io_registry.write(self, *args, **kwargs)


# Note the cls argument below; this is necessary here, but not in the Table
# reader because that reader assumes you are reading a Table or subclass.
# Here we have no idea what the subclass might be, so we state it explicitly.

# Well, try plan B: Use a closure to return a read_nddata_fits that returns
# the correct class.
def reader_for_class(cls):
    def read_nddata_fits(input, hdu=None):
        """
        Read FITS file into object that implements NDDataBase interface.

        Parameters
        ----------

        input : str or file-like object or compatible `astropy.io.fits` HDU object
            If a string, the filename to read the table from. If a file object, or
            a compatible HDU object, the object to extract the table from. The
            following `astropy.io.fits` HDU objects can be used as input:
            - :class:`~astropy.io.fits.hdu.table.TableHDU`
            - :class:`~astropy.io.fits.hdu.table.BinTableHDU`
            - :class:`~astropy.io.fits.hdu.table.GroupsHDU`
            - :class:`~astropy.io.fits.hdu.hdulist.HDUList`
        hdu : int or str, optional
            The HDU to read the table from.
        """
        def _is_nddata_hdu(hdu):
            return (isinstance(hdu, fits.ImageHDU) or
                    (isinstance(hdu, fits.PrimaryHDU) and
                     hdu.header['naxis'] > 0))
        if isinstance(input, fits.HDUList):
            nddata_hdus = OrderedDict()
            for idx, an_hdu in enumerate(input):
                if _is_nddata_hdu(an_hdu):
                    nddata_hdus[idx] = an_hdu
            if len(nddata_hdus) > 1:
                if hdu is None:
                    warnings.warn("hdu= was not specified but multiple tables"
                                  " are present, reading in first available"
                                  " table (hdu={0})".format(fits.util.first(nddata_hdus)),
                                  AstropyUserWarning)
                    hdu = fits.util.first(nddata_hdus)

                    # Grabbed the stuff below directly from table fits reader
                    # hdu might not be an integer, so we first need to convert it
                    # to the correct HDU index
                    hdu = input.index_of(hdu)

                    if hdu in nddata_hdus:
                        nddata = nddata_hdus[hdu]
                    else:
                        raise ValueError("No table found in hdu={0}".format(hdu))
            elif len(nddata_hdus) == 1:
                nddata = nddata_hdus[fits.util.first(nddata_hdus)]
            else:
                raise ValueError("No NDData-like extension found")

        elif _is_nddata_hdu(input):
            nddata = input
        else:
            # input is either a string or a file-like object, hopefully.
            hdulist = fits.open(input)
            try:
                return read_nddata_fits(hdulist, hdu=hdu)
            finally:
                hdulist.close()

        return cls(nddata.data, meta=nddata.header)
    return read_nddata_fits


def write_nddata_fits(input, output, overwrite=False):
    """
    Write an NDData-like object to a FITS file.

    Parameters
    ----------
    input : object that implements NDData interface
        The object whose data is to be written.
    output : str
        The filename to write the table to.
    overwrite : bool
        Whether to overwrite any existing file without warning.
    """

    # Thank you, io.fits.connect.write, for this snippet.
    # Check if output file already exists
    if isinstance(output, six.string_types) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    hdu = fits.PrimaryHDU(data=input.data, header=fits.Header(input.meta))
    if not isinstance(input.meta, fits.Header):
        # io.fits.Header doesn't use values when initializing from a dict,
        # it just creats the keys, so we'll enter the values manually.
        for k, v in six.iteritems(input.meta):
            hdu.header[k] = v

    hdu.writeto(output)
