# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os
import re
import warnings

import numpy as np

from ..registry.core import BaseIO
from ... import log
from ... import units as u
from ...extern.six import string_types
from ...table import Table
from ...nddata import NDData
from ...utils import OrderedDict
from ...utils.exceptions import AstropyUserWarning

from .hdu.base import _ValidHDU
from .hdu.table import _TableLikeHDU
from . import HDUList, Header
from .hdu.hdulist import fitsopen as fits_open
from .util import first


# FITS file signature as per RFC 4047
FITS_SIGNATURE = (b"\x53\x49\x4d\x50\x4c\x45\x20\x20\x3d\x20\x20\x20\x20\x20"
                  b"\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20"
                  b"\x20\x54")

# Keywords to remove for all tables that are read in
REMOVE_KEYWORDS = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'PCOUNT', 'GCOUNT', 'TFIELDS']

# Column-specific keywords
COLUMN_KEYWORDS = ['TFORM[0-9]+',
                   'TBCOL[0-9]+',
                   'TSCAL[0-9]+',
                   'TZERO[0-9]+',
                   'TNULL[0-9]+',
                   'TTYPE[0-9]+',
                   'TUNIT[0-9]+',
                   'TDISP[0-9]+',
                   'TDIM[0-9]+',
                   'THEAP']


def is_column_keyword(keyword):
    for c in COLUMN_KEYWORDS:
        if re.match(c, keyword) is not None:
            return True
    return False


def header_to_meta(header):

    meta = OrderedDict()

    for key, value, comment in header.cards:

        if key in ['COMMENT', 'HISTORY']:
            if key in meta:
                meta[key].append(value)
            else:
                meta[key] = [value]

        elif key in meta:  # key is duplicate

            if isinstance(meta[key], list):
                meta[key].append(value)
            else:
                meta[key] = [meta[key], value]

        elif (is_column_keyword(key.upper()) or
              key.upper() in REMOVE_KEYWORDS):

            pass

        else:

            meta[key] = value

    return meta


def meta_to_header(meta):

    header = Header()

    for key, value in meta.items():

        if is_column_keyword(key.upper()) or key.upper() in REMOVE_KEYWORDS:

            warnings.warn(
                "Meta-data keyword {0} will be ignored since it conflicts "
                "with a FITS reserved keyword".format(key), AstropyUserWarning)

        if isinstance(value, list):
            for item in value:
                try:
                    header.append((key, item))
                except ValueError:
                    warnings.warn(
                        "Attribute `{0}` of type {1} cannot be written to "
                        "FITS files - skipping".format(key, type(value)),
                        AstropyUserWarning)
        else:
            try:
                header[key] = value
            except ValueError:
                warnings.warn(
                    "Attribute `{0}` of type {1} cannot be written to FITS "
                    "files - skipping".format(key, type(value)),
                    AstropyUserWarning)

    return header


def _is_fits(origin, filepath, fileobj, *args, **kwargs):
    """
    Determine whether `origin` is a FITS file.

    Parameters
    ----------
    origin : str or readable file-like object
        Path or file object containing a potential FITS file.

    Returns
    -------
    is_fits : bool
        Returns `True` if the given file is a FITS file.
    """
    if fileobj is not None:
        pos = fileobj.tell()
        sig = fileobj.read(30)
        fileobj.seek(pos)
        return sig == FITS_SIGNATURE
    elif filepath is not None:
        if filepath.lower().endswith(('.fits', '.fits.gz', '.fit', '.fit.gz')):
            return True
    elif isinstance(args[0], (HDUList, _ValidHDU)):
        return True
    else:
        return False


class FITSTableIO(BaseIO):

    _format_name = 'fits'
    _supported_class = Table

    @staticmethod
    def identify(origin, filepath, fileobj, *args, **kwargs):
        return _is_fits(origin, filepath, fileobj, *args, **kwargs)

    @staticmethod
    def read(input, hdu=None):
        """
        Read a Table object from an FITS file

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

        if isinstance(input, HDUList):

            # Parse all table objects
            tables = OrderedDict()
            for ihdu, hdu_item in enumerate(input):
                if isinstance(hdu_item, (TableHDU, BinTableHDU, GroupsHDU)):
                    tables[ihdu] = hdu_item

            if len(tables) > 1:
                if hdu is None:
                    warnings.warn("hdu= was not specified but multiple tables"
                                  " are present, reading in first available"
                                  " table (hdu={0})".format(first(tables)),
                                  AstropyUserWarning)
                    hdu = first(tables)

                # hdu might not be an integer, so we first need to convert it
                # to the correct HDU index
                hdu = input.index_of(hdu)

                if hdu in tables:
                    table = tables[hdu]
                else:
                    raise ValueError("No table found in hdu={0}".format(hdu))

            elif len(tables) == 1:
                table = tables[first(tables)]
            else:
                raise ValueError("No table found")

        elif isinstance(input, _TableLikeHDU):

            table = input

        else:

            hdulist = fits_open(input)

            try:
                return Table.read(hdulist, hdu=hdu, format='fits')
            finally:
                hdulist.close()

        # Check if table is masked
        masked = False
        for col in table.columns:
            if col.null is not None:
                masked = True
                break

        # Convert to an astropy.table.Table object
        t = Table(table.data, masked=masked)

        # Copy over null values if needed
        if masked:
            for col in table.columns:
                if col.null is not None:
                    t[col.name].set_fill_value(col.null)
                    t[col.name].mask[t[col.name] == col.null] = True

        # Copy over units
        for col in table.columns:
            if col.unit is not None:
                t[col.name].unit = u.Unit(
                    col.unit, format='fits', parse_strict='warn')

        # TODO: deal properly with unsigned integers

        t.meta.update(header_to_meta(table.header))

        return t

    @staticmethod
    def write(input, output, overwrite=False):
        """
        Write a Table object to a FITS file

        Parameters
        ----------
        input : Table
            The table to write out.
        output : str
            The filename to write the table to.
        overwrite : bool
            Whether to overwrite any existing file without warning.
        """

        # Check if output file already exists
        if isinstance(output, string_types) and os.path.exists(output):
            if overwrite:
                os.remove(output)
            else:
                raise IOError("File exists: {0}".format(output))

        # Create a new HDU object
        if input.masked:
            table_hdu = BinTableHDU(np.array(input.filled()))
            for col in table_hdu.columns:
                # Binary FITS tables support TNULL *only* for integer data columns
                # TODO: Determine a schema for handling non-integer masked columns
                # in FITS (if at all possible)
                int_formats = ('B', 'I', 'J', 'K')
                if not (col.format in int_formats or
                        col.format.p_format in int_formats):
                    continue

                # The astype is necessary because if the string column is less
                # than one character, the fill value will be N/A by default which
                # is too long, and so no values will get masked.
                fill_value = input[col.name].get_fill_value()

                col.null = fill_value.astype(input[col.name].dtype)
        else:
            table_hdu = BinTableHDU(np.array(input))

        # Set units for output HDU
        for col in table_hdu.columns:
            if input[col.name].unit is not None:
                col.unit = input[col.name].unit.to_string(format='fits')

        table_hdu.header.update(meta_to_header(input.meta))

        # Write out file
        table_hdu.writeto(output)


class FITSNDDataIO(BaseIO):

    _format_name = 'fits'
    _supported_class = NDData

    @staticmethod
    def identify(origin, filepath, fileobj, *args, **kwargs):
        return _identify_fits(origin, filepath, fileobj, *args, **kwargs)

    @staticmethod
    def read(input, hdu=None):
        """
        Read an :class:`~astropy.nddata.nddata.NDData` object from a FITS file.

        Parameters
        ----------
        input : str or file-like object or compatible `astropy.io.fits` HDU object
            If a string, the filename to read the data from. If a file object, or
            a compatible HDU object, the object to extract the table from. The
            following `astropy.io.fits` HDU objects can be used as input:
            - :class:`~astropy.io.fits.hdu.table.PrimaryHDU`
            - :class:`~astropy.io.fits.hdu.table.ImageHDU`
            - :class:`~astropy.io.fits.hdu.hdulist.HDUList`
        hdu : int or str, optional
            The HDU to read the data from.
        """

        if isinstance(input, HDUList):

            # Parse all table objects
            images = OrderedDict()
            for ihdu, hdu_item in enumerate(input):
                if isinstance(hdu_item, (PrimaryHDU, ImageHDU)):
                    images[ihdu] = hdu_item

            if len(images) > 1:
                if hdu is None:
                    warnings.warn("hdu= was not specified but multiple n-dimensional images"
                                  " are present, reading in first available"
                                  " image (hdu={0})".format(first(images)),
                                  AstropyUserWarning)
                    hdu = first(images)

                # hdu might not be an integer, so we first need to convert it
                # to the correct HDU index
                hdu = input.index_of(hdu)

                if hdu in tables:
                    images = images[hdu]
                else:
                    raise ValueError("No image found in hdu={0}".format(hdu))

            elif len(images) == 1:
                image = images[first(images)]
            else:
                raise ValueError("No table found")

        elif isinstance(input, (PrimaryHDU, ImageHDU)):

            image = input

        else:

            hdulist = fits_open(input)

            try:
                return NDData.read(hdulist, hdu=hdu, format='fits')
            finally:
                hdulist.close()

        # Convert to an astropy.table.Table object
        d = NDData(image.data)
        if 'BUNIT' in image.header:
            d.unit = image.header['BUNIT']
        d.meta.update(header_to_meta(image.header))

        return d

    @staticmethod
    def write(input, output, overwrite=False):
        """
        Write an :class:`~astropy.nddata.nddata.NDData` object to a FITS file.

        Parameters
        ----------
        input : Table
            The table to write out.
        output : str
            The filename to write the table to.
        overwrite : bool
            Whether to overwrite any existing file without warning.
        """

        # Check if output file already exists
        if isinstance(output, string_types) and os.path.exists(output):
            if overwrite:
                os.remove(output)
            else:
                raise IOError("File exists: {0}".format(output))

        # Create a new HDU object
        image_hdu = PrimaryHDU(np.array(input))
        if input.unit is not None:
            image_hdu.header['BUNIT'] = input.unit.to_string(format='fits')
        image_hdu.header.update(meta_to_header(input.meta))
        image_hdu.writeto(output)
