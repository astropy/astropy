# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os
import re
import warnings
from collections import OrderedDict

import numpy as np

from .. import registry as io_registry
from ... import units as u
from ... import log
from ...extern import six
from ...extern.six import string_types
from ...table import Table
from ...utils.exceptions import AstropyUserWarning
from astropy.units.format.fits import UnitScaleError

from ...nddata import NDData, StdDevUncertainty, UnknownUncertainty, NDIOMixin
from ...wcs import WCS

from . import (HDUList, TableHDU, BinTableHDU, GroupsHDU, PrimaryHDU, ImageHDU,
               Header)
from . import FITS_rec
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


def is_fits(origin, filepath, fileobj, *args, **kwargs):
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
    elif isinstance(args[0], (HDUList, TableHDU, BinTableHDU, GroupsHDU)):
        return True
    else:
        return False


def read_table_fits(input, hdu=None):
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

    elif isinstance(input, (TableHDU, BinTableHDU, GroupsHDU)):

        table = input

    else:

        hdulist = fits_open(input)

        try:
            return read_table_fits(hdulist, hdu=hdu)
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

    for key, value, comment in table.header.cards:

        if key in ['COMMENT', 'HISTORY']:
            if key in t.meta:
                t.meta[key].append(value)
            else:
                t.meta[key] = [value]

        elif key in t.meta:  # key is duplicate

            if isinstance(t.meta[key], list):
                t.meta[key].append(value)
            else:
                t.meta[key] = [t.meta[key], value]

        elif (is_column_keyword(key.upper()) or
              key.upper() in REMOVE_KEYWORDS):

            pass

        else:

            t.meta[key] = value

    # TODO: implement masking

    return t


def write_table_fits(input, output, overwrite=False):
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

    # Tables with mixin columns are not supported
    if input.has_mixin_columns:
        mixin_names = [name for name, col in input.columns.items()
                       if not isinstance(col, input.ColumnClass)]
        raise ValueError('cannot write table with mixin column(s) {0} to FITS'
                         .format(mixin_names))

    # Check if output file already exists
    if isinstance(output, string_types) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    # Create a new HDU object
    if input.masked:
        #float column's default mask value needs to be Nan
        for column in six.itervalues(input.columns):
            fill_value = column.get_fill_value()
            if column.dtype.kind == 'f' and np.allclose(fill_value, 1e20):
                column.set_fill_value(np.nan)

        fits_rec = FITS_rec.from_columns(np.array(input.filled()))
        table_hdu = BinTableHDU(fits_rec)
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
        fits_rec = FITS_rec.from_columns(np.array(input.filled()))
        table_hdu = BinTableHDU(fits_rec)

    # Set units for output HDU
    for col in table_hdu.columns:
        unit = input[col.name].unit
        if unit is not None:
            try:
                col.unit = unit.to_string(format='fits')
            except UnitScaleError:
                scale = unit.scale
                raise UnitScaleError(
                    "The column '{0}' could not be stored in FITS format "
                    "because it has a scale '({1})' that "
                    "is not recognized by the FITS standard. Either scale "
                    "the data or change the units.".format(col.name, str(scale)))
            except ValueError:
                warnings.warn(
                    "The unit '{0}' could not be saved to FITS format".format(
                        unit.to_string()), AstropyUserWarning)

    for key, value in input.meta.items():

        if is_column_keyword(key.upper()) or key.upper() in REMOVE_KEYWORDS:

            warnings.warn(
                "Meta-data keyword {0} will be ignored since it conflicts "
                "with a FITS reserved keyword".format(key), AstropyUserWarning)

        if isinstance(value, list):
            for item in value:
                try:
                    table_hdu.header.append((key, item))
                except ValueError:
                    warnings.warn(
                        "Attribute `{0}` of type {1} cannot be written to "
                        "FITS files - skipping".format(key, type(value)),
                        AstropyUserWarning)
        else:
            try:
                table_hdu.header[key] = value
            except ValueError:
                warnings.warn(
                    "Attribute `{0}` of type {1} cannot be written to FITS "
                    "files - skipping".format(key, type(value)),
                    AstropyUserWarning)

    # Write out file
    table_hdu.writeto(output)

io_registry.register_reader('fits', Table, read_table_fits)
io_registry.register_writer('fits', Table, write_table_fits)
io_registry.register_identifier('fits', Table, is_fits)


def read_data_fits(filename, ext_data=0, ext_meta=0, ext_mask='mask',
                   ext_uncert='uncert', kw_unit='bunit', copy=False,
                   dtype=None, **kwargs_for_open):
    """
    Read data from a FITS file and wrap the contents in a \
    `~astropy.nddata.NDData`.

    Parameters
    ----------
    filename : str and other types
        see :func:`~astropy.io.fits.open` what possible types are allowed.

    ext_data, ext_meta, ext_mask, ext_uncert : str or int, optional
        Extensions from which to read ``data``, ``meta``, ``mask`` and
        ``uncertainty``.
        Default is ``0`` (data), ``0`` (meta), ``'mask'`` (mask) and
        ``'uncert'`` (uncertainty).

    kw_unit : str or None, optional
        The header keyword which translates to the unit for the data. Set it
        to ``None`` if parsing the unit results in a ValueError during reading.
        Default is ``'bunit'``.

    copy : bool, optional
        Copy the data while creating the `~astropy.nddata.NDData` instance?
        Default is ``False``.

    dtype : `numpy.dtype`-like or None, optional
        If not ``None`` the data array is converted to this dtype before
        returning. See `numpy.ndarray.astype` for more details.
        Default is ``None``.

    kwargs_for_open :
        Additional keyword arguments that are passed to
        :func:`~astropy.io.fits.open` (not all of them might be possible).
    """

    # Hardcoded values to get additional information about mask and uncertainty
    kw_hdr_masktype = 'boolean mask'
    kw_hdr_uncerttype = {
        'standard deviation uncertainty': StdDevUncertainty,
        'unknown uncertainty type': UnknownUncertainty}

    with fits_open(filename, mode='readonly', **kwargs_for_open) as hdus:
        # Read the data and meta from the specified extensions
        data = hdus[ext_data].data
        if dtype is not None:
            data = data.astype(dtype)
        meta = hdus[ext_meta].header

        # Read the mask and uncertainty from the specified extensions but
        # silently fail if the extension does not exist.
        mask = None
        if ext_mask in hdus:
            mask = hdus[ext_mask].data
            # Convert it to boolean array?
            if kw_hdr_masktype in hdus[ext_mask].header.get('comment', []):
                mask = mask.astype(bool)

        uncertainty = None
        if ext_uncert in hdus:
            uncertainty = hdus[ext_uncert].data

            hdr = hdus[ext_uncert].header
            # Get the required class for the uncertainty
            cls = (kw_hdr_uncerttype[kw] for kw in kw_hdr_uncerttype
                   if kw in hdr.get('comment', []))
            cls = next(cls, UnknownUncertainty)

            # Get the unit for the uncertainty if present
            unit_ = hdr[kw_unit].lower() if kw_unit in hdr else None

            # Don't copy it here, if a copy is required do it when creating
            # NDData.
            uncertainty = cls(uncertainty, unit=unit_, copy=False)

        # Load unit and wcs from header
        unit = None
        if kw_unit is not None and kw_unit in meta:
            try:
                unit = u.Unit(meta[kw_unit])
            except ValueError:
                # ValueError is raised if the unit isn't convertible to an
                # astropy unit. Maybe they tried all-uppercase: maybe lowercase
                # will work
                unit = u.Unit(meta[kw_unit].lower())
        wcs = WCS(meta)

    # Just create an NDData instance: This will be upcast to the appropriate
    # class
    return NDData(data, meta=meta, mask=mask, uncertainty=uncertainty,
                  wcs=wcs, unit=unit, copy=copy)


def write_data_fits(ndd, filename, ext_mask='mask', ext_uncert='uncert',
                    kw_unit='bunit', **kwargs_for_write):
    """
    Take an `~astropy.nddata.NDData`-like object and save it as FITS file.

    Parameters
    ----------
    ndd : `astropy.nddata.NDData`-like
        The data which is to be saved. Must not be given when this function
        is called through the ``NDData.write``-method!

    filename : str
        The filename for the newly written file.

    ext_mask, ext_uncert : str or int, optional
        Extensions to which ``mask`` and ``uncertainty`` are written.
        Default is ``'mask'`` (mask) and ``'uncert'`` (uncertainty).

    kwargs_for_write :
        Additional keyword arguments that are passed to
        :func:`~astropy.io.fits.HDUList.writeto` (not all of them might be
        possible).

    Notes
    -----
    The ``data`` and ``meta`` are always written to the PrimaryHDU (extension
    number ``0``).
    """
    # Comment card strings to allow roundtripping (must be identical to read!)
    kw_hdr_masktype = 'boolmask'
    kw_hdr_uncerttype = {
        StdDevUncertainty: 'standard deviation uncertainty',
        UnknownUncertainty: 'unknown uncertainty type'}

    # Copy or convert the meta to a FITS header
    if isinstance(ndd.meta, Header):
        header = ndd.meta.copy()
    else:
        header = Header(ndd.meta.items())

    # Update the (copied) header (unit, wcs)
    if ndd.unit is not None:
        header[kw_unit] = ndd.unit.to_string()
    elif kw_unit in header:
        del header[kw_unit]

    if ndd.wcs is not None:
        try:
            header.update(ndd.wcs.to_header())
        except AttributeError:
            # wcs has no to_header method
            # FIXME: Implement this if other wcs objects should be allowed.
            log.info("the wcs cannot be converted to header information.")

    # Create a HDUList containing data
    hdus = [PrimaryHDU(ndd.data, header=header)]

    # And append mask to the HDUList (if present)
    try:
        # Convert mask to uint8 and set a keyword so that the opener knows
        # that it was a boolean mask and can convert it back again.
        if ndd.mask.dtype == 'bool':
            hdr = Header()
            hdr.add_comment(kw_hdr_masktype)
            hdus.append(ImageHDU(ndd.mask.astype(np.uint8), header=hdr,
                                 name=ext_mask))
        else:
            hdus.append(ImageHDU(ndd.mask, name=ext_mask))
    except AttributeError:
        # Either no mask or mask had no dtype
        pass

    # And append the uncertainty (if present)
    try:
        # We need to save the uncertainty_type and the unit of the uncertainty
        # so that the uncertainty can be completly recovered.
        hdr = Header()

        # Save the class of the uncertainty
        if ndd.uncertainty.__class__ in kw_hdr_uncerttype:
            hdr.add_comment(kw_hdr_uncerttype[ndd.uncertainty.__class__])

        # Save the unit of the uncertainty if it differs from the nddata
        # TODO: This comparison only works correctly for StdDevUncertainty...
        if ndd.uncertainty.unit != ndd.unit:
            hdr[kw_unit] = ndd.uncertainty.unit.to_string()

        hdus.append(ImageHDU(ndd.uncertainty.array, header=hdr,
                             name=ext_uncert))
    except AttributeError:
        # Either no uncertainty or no uncertainty array, unit or
        # uncertainty_type. Should not be possible because everything that
        # doesn't look like an NDUUncertainty is converted to one.
        pass

    # Convert to HDUList and write it to the file.
    with HDUList(hdus) as hdulist:
        hdulist.writeto(filename, **kwargs_for_write)


# TODO: Register reader and writer WITHOUT identifier (for now...)
io_registry.register_reader('simple_fits', NDIOMixin, read_data_fits)
io_registry.register_writer('simple_fits', NDIOMixin, write_data_fits)
