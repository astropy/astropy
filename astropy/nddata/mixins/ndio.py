# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the I/O mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ...io import registry as io_registry
from ...io import fits
from ...wcs import WCS
from .. import NDData, StdDevUncertainty, UnknownUncertainty

__all__ = ['NDIOMixin']


class NDIOMixin(object):
    """
    Mixin class to connect NDData to the astropy input/output registry.

    This mixin adds two methods to its subclasses, ``read`` and ``write``.
    """

    @classmethod
    def read(cls, *args, **kwargs):
        """
        Read and parse gridded N-dimensional data and return as an
        NDData-derived object.

        This function provides the NDDataBase interface to the astropy unified
        I/O layer.  This allows easily reading a file in the supported data
        formats.
        """
        return io_registry.read(cls, *args, **kwargs)

    def write(self, *args, **kwargs):
        """
        Write a gridded N-dimensional data object out in specified format.

        This function provides the NDDataBase interface to the astropy unified
        I/O layer.  This allows easily writing a file in the supported data
        formats.
        """
        io_registry.write(self, *args, **kwargs)


def read_from_fits(filename):
    # Variables
    ext_data = 0
    ext_meta = 0
    ext_mask = 'mask'
    ext_uncert = 'uncert'
    kw_unit = 'bunit'
    kw_mask_type = 'boolmask'
    kw_uncert_type = 'stddev'

    with fits.open(filename) as hdus:
        # Read the data from the primary hdu
        data = hdus[ext_data].data
        # Meta is also taken from the primary
        meta = hdus[ext_meta].header
        # Mask is rad from the MASK extension ImageHDU
        mask = None
        if ext_mask in hdus:
            mask = hdus[ext_mask].data
            # Check if it's a boolean mask ('boolmask' is set in the extension
            # header) and convert it accordingly
            if kw_mask_type in hdus[ext_mask].header:
                mask = mask.astype(bool)

        # Same for the uncertainty
        uncertainty = None
        if ext_uncert in hdus:
            uncertainty = hdus[ext_uncert].data
            # Now lets check if its standard deviation and convert it
            # accordingly
            if kw_uncert_type in hdus[ext_uncert].header:
                cls = StdDevUncertainty
            else:
                cls = UnknownUncertainty

            if kw_unit in hdus[ext_uncert].header:
                uncert_unit = hdus[ext_uncert].header[kw_unit].lower()
            else:
                uncert_unit = None
            uncertainty = cls(uncertainty, unit=uncert_unit)

        # Check if a unit is given. Fits standard say 'BUNIT'-keyword
        # Convert to lowercase so that astropy.units.Unit knows which it is.
        unit = str(meta[kw_unit]).lower() if kw_unit in meta else None

        # If it had a valid primary header it will be convertable to an
        # astropy.wcs.WCS object:
        wcs = WCS(meta)

    # Closing the file and returning a NDData instance. This will be upcast if
    # this function was called by a subclass of NDIOMixin.read
    return NDData(data, meta=meta, mask=mask, uncertainty=uncertainty,
                  wcs=wcs, unit=unit)


def write_to_fits(ndd, filename):
    # Variables
    # ext_data = 0  # data can only be written into primary
    # ext_meta = 0  # header too
    ext_mask = 'mask'
    ext_uncert = 'uncert'
    kw_unit = 'bunit'
    kw_mask_type = 'boolmask'
    kw_uncert_type = 'stddev'

    # Convert to a fits.Header or copy (if appropriate) the meta.
    if isinstance(ndd.meta, fits.Header):
        header = ndd.meta.copy()
    else:
        header = fits.Header(ndd.meta.items())

    # Update any WCS changes to header and set/reset the unit
    if ndd.unit is not None:
        header[kw_unit] = ndd.unit.to_string()
    elif kw_unit in header:
        # Delete the unit from the header if the nddata had no unit
        del header[kw_unit]

    if ndd.wcs is not None:
        header.update(ndd.wcs.to_header())

    # Create a HDUList containing data
    hdus = [fits.PrimaryHDU(ndd.data, header=header)]

    # And append mask and uncertainty to the HDUList (if present)
    try:
        # I haven't figured out how to write boolean arrays to FITS so
        # I convert it to uint8 and set a keyword so that the opener knows
        # that it was a boolean mask and can convert it to one again.
        if ndd.mask.dtype == 'bool':
            hdr = fits.Header([(kw_mask_type, 'True')])
            hdus.append(fits.ImageHDU(ndd.mask.astype(np.uint8),
                                      header=hdr,
                                      name=ext_mask))
        else:
            # It had a dtype so just save it as ImageHDU.
            hdus.append(fits.ImageHDU(ndd.mask,
                                      name=ext_mask))
    except AttributeError:
        # Either no mask or mask had no dtype
        pass

    # Same for uncertainty but without boolean shortcut we would not
    # expect uncertainty to contain boolean data
    try:
        # We need to save the uncertainty_type and the unit of the uncertainty
        # so that the uncertainty can be recovered.
        if ndd.uncertainty.uncertainty_type == 'std':
            hdr = fits.Header([(kw_uncert_type, 'True')])
        else:
            hdr = fits.Header()
        if ndd.uncertainty.unit != ndd.unit:
            hdr[kw_unit] = ndd.uncertainty.unit.to_string()
        hdus.append(fits.ImageHDU(ndd.uncertainty.array,
                                  header=hdr,
                                  name=ext_uncert))
    except AttributeError:
        # Either no uncertainty or no uncertainty array, unit or
        # uncertainty_type
        pass

    # TODO: Maybe add another except in case the uncertainty array wasn't a
    # numpy array or something that could be saved in an ImageHDU. But I'll
    # assume for now that letting the exception rise is better for now.

    # Converting to HDUList and writing.
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(filename)

    # do I have to do something with hdulist here to avoid open file handles?


# Register reader and writer WITHOUT identifier (for now)
io_registry.register_reader('simple_fits', NDIOMixin, read_from_fits)
io_registry.register_writer('simple_fits', NDIOMixin, write_to_fits)
