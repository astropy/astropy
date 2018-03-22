# Licensed under a 3-clause BSD style license - see PYFITS.rst

"""
Convenience functions
=====================

The functions in this module provide shortcuts for some of the most basic
operations on FITS files, such as reading and updating the header.  They are
included directly in the 'astropy.io.fits' namespace so that they can be used
like::

    astropy.io.fits.getheader(...)

These functions are primarily for convenience when working with FITS files in
the command-line interpreter.  If performing several operations on the same
file, such as in a script, it is better to *not* use these functions, as each
one must open and re-parse the file.  In such cases it is better to use
:func:`astropy.io.fits.open` and work directly with the
:class:`astropy.io.fits.HDUList` object and underlying HDU objects.

Several of the convenience functions, such as `getheader` and `getdata` support
special arguments for selecting which extension HDU to use when working with a
multi-extension FITS file.  There are a few supported argument formats for
selecting the extension.  See the documentation for `getdata` for an
explanation of all the different formats.

.. warning::
    All arguments to convenience functions other than the filename that are
    *not* for selecting the extension HDU should be passed in as keyword
    arguments.  This is to avoid ambiguity and conflicts with the
    extension arguments.  For example, to set NAXIS=1 on the Primary HDU:

    Wrong::

        astropy.io.fits.setval('myimage.fits', 'NAXIS', 1)

    The above example will try to set the NAXIS value on the first extension
    HDU to blank.  That is, the argument '1' is assumed to specify an extension
    HDU.

    Right::

        astropy.io.fits.setval('myimage.fits', 'NAXIS', value=1)

    This will set the NAXIS keyword to 1 on the primary HDU (the default).  To
    specify the first extension HDU use::

        astropy.io.fits.setval('myimage.fits', 'NAXIS', value=1, ext=1)

    This complexity arises out of the attempt to simultaneously support
    multiple argument formats that were used in past versions of PyFITS.
    Unfortunately, it is not possible to support all formats without
    introducing some ambiguity.  A future Astropy release may standardize
    around a single format and officially deprecate the other formats.
"""


import operator
import os
import warnings

import numpy as np

from .diff import FITSDiff, HDUDiff
from .file import FILE_MODES, _File
from .hdu.base import _BaseHDU, _ValidHDU
from .hdu.hdulist import fitsopen, HDUList
from .hdu.image import PrimaryHDU, ImageHDU
from .hdu.table import BinTableHDU
from .header import Header
from .util import fileobj_closed, fileobj_name, fileobj_mode, _is_int
from ...units import Unit
from ...units.format.fits import UnitScaleError
from ...units import Quantity
from ...utils.exceptions import AstropyUserWarning
from ...utils.decorators import deprecated_renamed_argument


__all__ = ['getheader', 'getdata', 'getval', 'setval', 'delval', 'writeto',
           'append', 'update', 'info', 'tabledump', 'tableload',
           'table_to_hdu', 'printdiff']


def getheader(filename, *args, **kwargs):
    """
    Get the header from an extension of a FITS file.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to get header from.  If an opened file object, its mode
        must be one of the following rb, rb+, or ab+).

    ext, extname, extver
        The rest of the arguments are for extension specification.  See the
        `getdata` documentation for explanations/examples.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.

    Returns
    -------
    header : `Header` object
    """

    mode, closed = _get_file_mode(filename)
    hdulist, extidx = _getext(filename, mode, *args, **kwargs)
    try:
        hdu = hdulist[extidx]
        header = hdu.header
    finally:
        hdulist.close(closed=closed)

    return header


def getdata(filename, *args, header=None, lower=None, upper=None, view=None,
            **kwargs):
    """
    Get the data from an extension of a FITS file (and optionally the
    header).

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to get data from.  If opened, mode must be one of the
        following rb, rb+, or ab+.

    ext
        The rest of the arguments are for extension specification.
        They are flexible and are best illustrated by examples.

        No extra arguments implies the primary header::

            getdata('in.fits')

        By extension number::

            getdata('in.fits', 0)      # the primary header
            getdata('in.fits', 2)      # the second extension
            getdata('in.fits', ext=2)  # the second extension

        By name, i.e., ``EXTNAME`` value (if unique)::

            getdata('in.fits', 'sci')
            getdata('in.fits', extname='sci')  # equivalent

        Note ``EXTNAME`` values are not case sensitive

        By combination of ``EXTNAME`` and EXTVER`` as separate
        arguments or as a tuple::

            getdata('in.fits', 'sci', 2)  # EXTNAME='SCI' & EXTVER=2
            getdata('in.fits', extname='sci', extver=2)  # equivalent
            getdata('in.fits', ('sci', 2))  # equivalent

        Ambiguous or conflicting specifications will raise an exception::

            getdata('in.fits', ext=('sci',1), extname='err', extver=2)

    header : bool, optional
        If `True`, return the data and the header of the specified HDU as a
        tuple.

    lower, upper : bool, optional
        If ``lower`` or ``upper`` are `True`, the field names in the
        returned data object will be converted to lower or upper case,
        respectively.

    view : ndarray, optional
        When given, the data will be returned wrapped in the given ndarray
        subclass by calling::

           data.view(view)

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.

    Returns
    -------
    array : array, record array or groups data object
        Type depends on the type of the extension being referenced.

        If the optional keyword ``header`` is set to `True`, this
        function will return a (``data``, ``header``) tuple.
    """

    mode, closed = _get_file_mode(filename)

    hdulist, extidx = _getext(filename, mode, *args, **kwargs)
    try:
        hdu = hdulist[extidx]
        data = hdu.data
        if data is None and extidx == 0:
            try:
                hdu = hdulist[1]
                data = hdu.data
            except IndexError:
                raise IndexError('No data in this HDU.')
        if data is None:
            raise IndexError('No data in this HDU.')
        if header:
            hdr = hdu.header
    finally:
        hdulist.close(closed=closed)

    # Change case of names if requested
    trans = None
    if lower:
        trans = operator.methodcaller('lower')
    elif upper:
        trans = operator.methodcaller('upper')
    if trans:
        if data.dtype.names is None:
            # this data does not have fields
            return
        if data.dtype.descr[0][0] == '':
            # this data does not have fields
            return
        data.dtype.names = [trans(n) for n in data.dtype.names]

    # allow different views into the underlying ndarray.  Keep the original
    # view just in case there is a problem
    if isinstance(view, type) and issubclass(view, np.ndarray):
        data = data.view(view)

    if header:
        return data, hdr
    else:
        return data


def getval(filename, keyword, *args, **kwargs):
    """
    Get a keyword's value from a header in a FITS file.

    Parameters
    ----------
    filename : file path, file object, or file like object
        Name of the FITS file, or file object (if opened, mode must be
        one of the following rb, rb+, or ab+).

    keyword : str
        Keyword name

    ext, extname, extver
        The rest of the arguments are for extension specification.
        See `getdata` for explanations/examples.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
        *Note:* This function automatically specifies ``do_not_scale_image_data
        = True`` when opening the file so that values can be retrieved from the
        unmodified header.

    Returns
    -------
    keyword value : str, int, or float
    """

    if 'do_not_scale_image_data' not in kwargs:
        kwargs['do_not_scale_image_data'] = True

    hdr = getheader(filename, *args, **kwargs)
    return hdr[keyword]


def setval(filename, keyword, *args, value=None, comment=None, before=None,
           after=None, savecomment=False, **kwargs):
    """
    Set a keyword's value from a header in a FITS file.

    If the keyword already exists, it's value/comment will be updated.
    If it does not exist, a new card will be created and it will be
    placed before or after the specified location.  If no ``before`` or
    ``after`` is specified, it will be appended at the end.

    When updating more than one keyword in a file, this convenience
    function is a much less efficient approach compared with opening
    the file for update, modifying the header, and closing the file.

    Parameters
    ----------
    filename : file path, file object, or file like object
        Name of the FITS file, or file object If opened, mode must be update
        (rb+).  An opened file object or `~gzip.GzipFile` object will be closed
        upon return.

    keyword : str
        Keyword name

    value : str, int, float, optional
        Keyword value (default: `None`, meaning don't modify)

    comment : str, optional
        Keyword comment, (default: `None`, meaning don't modify)

    before : str, int, optional
        Name of the keyword, or index of the card before which the new card
        will be placed.  The argument ``before`` takes precedence over
        ``after`` if both are specified (default: `None`).

    after : str, int, optional
        Name of the keyword, or index of the card after which the new card will
        be placed. (default: `None`).

    savecomment : bool, optional
        When `True`, preserve the current comment for an existing keyword.  The
        argument ``savecomment`` takes precedence over ``comment`` if both
        specified.  If ``comment`` is not specified then the current comment
        will automatically be preserved  (default: `False`).

    ext, extname, extver
        The rest of the arguments are for extension specification.
        See `getdata` for explanations/examples.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
        *Note:* This function automatically specifies ``do_not_scale_image_data
        = True`` when opening the file so that values can be retrieved from the
        unmodified header.
    """

    if 'do_not_scale_image_data' not in kwargs:
        kwargs['do_not_scale_image_data'] = True

    closed = fileobj_closed(filename)
    hdulist, extidx = _getext(filename, 'update', *args, **kwargs)
    try:
        if keyword in hdulist[extidx].header and savecomment:
            comment = None
        hdulist[extidx].header.set(keyword, value, comment, before, after)
    finally:
        hdulist.close(closed=closed)


def delval(filename, keyword, *args, **kwargs):
    """
    Delete all instances of keyword from a header in a FITS file.

    Parameters
    ----------

    filename : file path, file object, or file like object
        Name of the FITS file, or file object If opened, mode must be update
        (rb+).  An opened file object or `~gzip.GzipFile` object will be closed
        upon return.

    keyword : str, int
        Keyword name or index

    ext, extname, extver
        The rest of the arguments are for extension specification.
        See `getdata` for explanations/examples.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
        *Note:* This function automatically specifies ``do_not_scale_image_data
        = True`` when opening the file so that values can be retrieved from the
        unmodified header.
    """

    if 'do_not_scale_image_data' not in kwargs:
        kwargs['do_not_scale_image_data'] = True

    closed = fileobj_closed(filename)
    hdulist, extidx = _getext(filename, 'update', *args, **kwargs)
    try:
        del hdulist[extidx].header[keyword]
    finally:
        hdulist.close(closed=closed)


@deprecated_renamed_argument('clobber', 'overwrite', '2.0')
def writeto(filename, data, header=None, output_verify='exception',
            overwrite=False, checksum=False):
    """
    Create a new FITS file using the supplied data/header.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened in a writeable binary
        mode such as 'wb' or 'ab+'.

    data : array, record array, or groups data object
        data to write to the new file

    header : `Header` object, optional
        the header associated with ``data``. If `None`, a header
        of the appropriate type is created for the supplied data. This
        argument is optional.

    output_verify : str
        Output verification option.  Must be one of ``"fix"``, ``"silentfix"``,
        ``"ignore"``, ``"warn"``, or ``"exception"``.  May also be any
        combination of ``"fix"`` or ``"silentfix"`` with ``"+ignore"``,
        ``+warn``, or ``+exception" (e.g. ``"fix+warn"``).  See :ref:`verify`
        for more info.

    overwrite : bool, optional
        If ``True``, overwrite the output file if it exists. Raises an
        ``OSError`` if ``False`` and the output file exists. Default is
        ``False``.

        .. versionchanged:: 1.3
           ``overwrite`` replaces the deprecated ``clobber`` argument.

    checksum : bool, optional
        If `True`, adds both ``DATASUM`` and ``CHECKSUM`` cards to the
        headers of all HDU's written to the file.
    """

    hdu = _makehdu(data, header)
    if hdu.is_image and not isinstance(hdu, PrimaryHDU):
        hdu = PrimaryHDU(data, header=header)
    hdu.writeto(filename, overwrite=overwrite, output_verify=output_verify,
                checksum=checksum)


def table_to_hdu(table, character_as_bytes=False):
    """
    Convert an `~astropy.table.Table` object to a FITS
    `~astropy.io.fits.BinTableHDU`.

    Parameters
    ----------
    table : astropy.table.Table
        The table to convert.
    character_as_bytes : bool
        Whether to return bytes for string columns when accessed from the HDU.
        By default this is `False` and (unicode) strings are returned, but for
        large tables this may use up a lot of memory.

    Returns
    -------
    table_hdu : `~astropy.io.fits.BinTableHDU`
        The FITS binary table HDU.
    """
    # Avoid circular imports
    from .connect import is_column_keyword, REMOVE_KEYWORDS

    # Header to store Time related metadata
    hdr = None

    # Not all tables with mixin columns are supported
    if table.has_mixin_columns:
        # Import is done here, in order to avoid it at build time as erfa is not
        # yet available then.
        from ...table.column import BaseColumn, Column
        from ...time import Time
        from .fitstime import time_to_fits

        # Only those columns which are instances of BaseColumn, Quantity or Time can
        # be written
        unsupported_cols = table.columns.not_isinstance((BaseColumn, Quantity, Time))
        if unsupported_cols:
            unsupported_names = [col.info.name for col in unsupported_cols]
            raise ValueError('cannot write table with mixin column(s) {0}'
                         .format(unsupported_names))

        time_cols = table.columns.isinstance(Time)
        if time_cols:
            table, hdr = time_to_fits(table)

    # Create a new HDU object
    if table.masked:
        # float column's default mask value needs to be Nan
        for column in table.columns.values():
            fill_value = column.get_fill_value()
            if column.dtype.kind == 'f' and np.allclose(fill_value, 1e20):
                column.set_fill_value(np.nan)

        # TODO: it might be better to construct the FITS table directly from
        # the Table columns, rather than go via a structured array.
        table_hdu = BinTableHDU.from_columns(np.array(table.filled()), header=hdr, character_as_bytes=True)
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
            fill_value = table[col.name].get_fill_value()

            col.null = fill_value.astype(table[col.name].dtype)
    else:
        table_hdu = BinTableHDU.from_columns(np.array(table.filled()), header=hdr, character_as_bytes=character_as_bytes)

    # Set units for output HDU
    for col in table_hdu.columns:
        unit = table[col.name].unit
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

            # Try creating a Unit to issue a warning if the unit is not FITS compliant
            Unit(col.unit, format='fits', parse_strict='warn')

    # Column-specific override keywords for coordinate columns
    coord_meta = table.meta.pop('__coordinate_columns__', {})
    for col_name, col_info in coord_meta.items():
        col = table_hdu.columns[col_name]
        # Set the column coordinate attributes from data saved earlier.
        # Note: have to set all three, even if we have no data.
        for attr in 'coord_type', 'coord_unit', 'time_ref_pos':
            setattr(col, attr, col_info.get(attr, None))

    for key, value in table.meta.items():
        if is_column_keyword(key.upper()) or key.upper() in REMOVE_KEYWORDS:
            warnings.warn(
                "Meta-data keyword {0} will be ignored since it conflicts "
                "with a FITS reserved keyword".format(key), AstropyUserWarning)

        # Convert to FITS format
        if key == 'comments':
            key = 'comment'

        if isinstance(value, list):
            for item in value:
                try:
                    table_hdu.header.append((key, item))
                except ValueError:
                    warnings.warn(
                        "Attribute `{0}` of type {1} cannot be added to "
                        "FITS Header - skipping".format(key, type(value)),
                        AstropyUserWarning)
        else:
            try:
                table_hdu.header[key] = value
            except ValueError:
                warnings.warn(
                    "Attribute `{0}` of type {1} cannot be added to FITS "
                    "Header - skipping".format(key, type(value)),
                    AstropyUserWarning)
    return table_hdu


def append(filename, data, header=None, checksum=False, verify=True, **kwargs):
    """
    Append the header/data to FITS file if filename exists, create if not.

    If only ``data`` is supplied, a minimal header is created.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened for update (rb+) unless it
        is a new file, then it must be opened for append (ab+).  A file or
        `~gzip.GzipFile` object opened for update will be closed after return.

    data : array, table, or group data object
        the new data used for appending

    header : `Header` object, optional
        The header associated with ``data``.  If `None`, an appropriate header
        will be created for the data object supplied.

    checksum : bool, optional
        When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards to the header
        of the HDU when written to the file.

    verify : bool, optional
        When `True`, the existing FITS file will be read in to verify it for
        correctness before appending.  When `False`, content is simply appended
        to the end of the file.  Setting ``verify`` to `False` can be much
        faster.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
    """

    name, closed, noexist_or_empty = _stat_filename_or_fileobj(filename)

    if noexist_or_empty:
        #
        # The input file or file like object either doesn't exits or is
        # empty.  Use the writeto convenience function to write the
        # output to the empty object.
        #
        writeto(filename, data, header, checksum=checksum, **kwargs)
    else:
        hdu = _makehdu(data, header)

        if isinstance(hdu, PrimaryHDU):
            hdu = ImageHDU(data, header)

        if verify or not closed:
            f = fitsopen(filename, mode='append')
            try:
                f.append(hdu)

                # Set a flag in the HDU so that only this HDU gets a checksum
                # when writing the file.
                hdu._output_checksum = checksum
            finally:
                f.close(closed=closed)
        else:
            f = _File(filename, mode='append')
            try:
                hdu._output_checksum = checksum
                hdu._writeto(f)
            finally:
                f.close()


def update(filename, data, *args, **kwargs):
    """
    Update the specified extension with the input data/header.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to update.  If opened, mode must be update (rb+).  An opened file
        object or `~gzip.GzipFile` object will be closed upon return.

    data : array, table, or group data object
        the new data used for updating

    header : `Header` object, optional
        The header associated with ``data``.  If `None`, an appropriate header
        will be created for the data object supplied.

    ext, extname, extver
        The rest of the arguments are flexible: the 3rd argument can be the
        header associated with the data.  If the 3rd argument is not a
        `Header`, it (and other positional arguments) are assumed to be the
        extension specification(s).  Header and extension specs can also be
        keyword arguments.  For example::

            update(file, dat, hdr, 'sci')  # update the 'sci' extension
            update(file, dat, 3)  # update the 3rd extension
            update(file, dat, hdr, 3)  # update the 3rd extension
            update(file, dat, 'sci', 2)  # update the 2nd SCI extension
            update(file, dat, 3, header=hdr)  # update the 3rd extension
            update(file, dat, header=hdr, ext=5)  # update the 5th extension

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
    """

    # The arguments to this function are a bit trickier to deal with than others
    # in this module, since the documentation has promised that the header
    # argument can be an optional positional argument.
    if args and isinstance(args[0], Header):
        header = args[0]
        args = args[1:]
    else:
        header = None
    # The header can also be a keyword argument--if both are provided the
    # keyword takes precedence
    header = kwargs.pop('header', header)

    new_hdu = _makehdu(data, header)

    closed = fileobj_closed(filename)

    hdulist, _ext = _getext(filename, 'update', *args, **kwargs)
    try:
        hdulist[_ext] = new_hdu
    finally:
        hdulist.close(closed=closed)


def info(filename, output=None, **kwargs):
    """
    Print the summary information on a FITS file.

    This includes the name, type, length of header, data shape and type
    for each extension.

    Parameters
    ----------
    filename : file path, file object, or file like object
        FITS file to obtain info from.  If opened, mode must be one of
        the following: rb, rb+, or ab+ (i.e. the file must be readable).

    output : file, bool, optional
        A file-like object to write the output to.  If ``False``, does not
        output to a file and instead returns a list of tuples representing the
        HDU info.  Writes to ``sys.stdout`` by default.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
        *Note:* This function sets ``ignore_missing_end=True`` by default.
    """

    mode, closed = _get_file_mode(filename, default='readonly')
    # Set the default value for the ignore_missing_end parameter
    if 'ignore_missing_end' not in kwargs:
        kwargs['ignore_missing_end'] = True

    f = fitsopen(filename, mode=mode, **kwargs)
    try:
        ret = f.info(output=output)
    finally:
        if closed:
            f.close()

    return ret


def printdiff(inputa, inputb, *args, **kwargs):
    """
    Compare two parts of a FITS file, including entire FITS files,
    FITS `HDUList` objects and FITS ``HDU`` objects.

    Parameters
    ----------
    inputa : str, `HDUList` object, or ``HDU`` object
        The filename of a FITS file, `HDUList`, or ``HDU``
        object to compare to ``inputb``.

    inputb : str, `HDUList` object, or ``HDU`` object
        The filename of a FITS file, `HDUList`, or ``HDU``
        object to compare to ``inputa``.

    ext, extname, extver
        Additional positional arguments are for extension specification if your
        inputs are string filenames (will not work if
        ``inputa`` and ``inputb`` are ``HDU`` objects or `HDUList` objects).
        They are flexible and are best illustrated by examples.  In addition
        to using these arguments positionally you can directly call the
        keyword parameters ``ext``, ``extname``.

        By extension number::

            printdiff('inA.fits', 'inB.fits', 0)      # the primary HDU
            printdiff('inA.fits', 'inB.fits', 2)      # the second extension
            printdiff('inA.fits', 'inB.fits', ext=2)  # the second extension

        By name, i.e., ``EXTNAME`` value (if unique). ``EXTNAME`` values are
        not case sensitive:

            printdiff('inA.fits', 'inB.fits', 'sci')
            printdiff('inA.fits', 'inB.fits', extname='sci')  # equivalent

        By combination of ``EXTNAME`` and ``EXTVER`` as separate
        arguments or as a tuple::

            printdiff('inA.fits', 'inB.fits', 'sci', 2)    # EXTNAME='SCI'
                                                           # & EXTVER=2
            printdiff('inA.fits', 'inB.fits', extname='sci', extver=2)
                                                           # equivalent
            printdiff('inA.fits', 'inB.fits', ('sci', 2))  # equivalent

        Ambiguous or conflicting specifications will raise an exception::

            printdiff('inA.fits', 'inB.fits',
                      ext=('sci', 1), extname='err', extver=2)

    kwargs
        Any additional keyword arguments to be passed to
        `~astropy.io.fits.FITSDiff`.

    Notes
    -----
    The primary use for the `printdiff` function is to allow quick print out
    of a FITS difference report and will write to ``sys.stdout``.
    To save the diff report to a file please use `~astropy.io.fits.FITSDiff`
    directly.
    """

    # Pop extension keywords
    extension = {key: kwargs.pop(key) for key in ['ext', 'extname', 'extver']
                 if key in kwargs}
    has_extensions = args or extension

    if isinstance(inputa, str) and has_extensions:
        # Use handy _getext to interpret any ext keywords, but
        # will need to close a if  fails
        modea, closeda = _get_file_mode(inputa)
        modeb, closedb = _get_file_mode(inputb)

        hdulista, extidxa = _getext(inputa, modea, *args, **extension)
        # Have to close a if b doesn't make it
        try:
            hdulistb, extidxb = _getext(inputb, modeb, *args, **extension)
        except Exception:
            hdulista.close(closed=closeda)
            raise

        try:
            hdua = hdulista[extidxa]
            hdub = hdulistb[extidxb]
            # See below print for note
            print(HDUDiff(hdua, hdub, **kwargs).report())

        finally:
            hdulista.close(closed=closeda)
            hdulistb.close(closed=closedb)

    # If input is not a string, can feed HDU objects or HDUList directly,
    # but can't currently handle extensions
    elif isinstance(inputa, _ValidHDU) and has_extensions:
        raise ValueError("Cannot use extension keywords when providing an "
                         "HDU object.")

    elif isinstance(inputa, _ValidHDU) and not has_extensions:
        print(HDUDiff(inputa, inputb, **kwargs).report())

    elif isinstance(inputa, HDUList) and has_extensions:
        raise NotImplementedError("Extension specification with HDUList "
                                  "objects not implemented.")

    # This function is EXCLUSIVELY for printing the diff report to screen
    # in a one-liner call, hence the use of print instead of logging
    else:
        print(FITSDiff(inputa, inputb, **kwargs).report())


@deprecated_renamed_argument('clobber', 'overwrite', '2.0')
def tabledump(filename, datafile=None, cdfile=None, hfile=None, ext=1,
              overwrite=False):
    """
    Dump a table HDU to a file in ASCII format.  The table may be
    dumped in three separate files, one containing column definitions,
    one containing header parameters, and one for table data.

    Parameters
    ----------
    filename : file path, file object or file-like object
        Input fits file.

    datafile : file path, file object or file-like object, optional
        Output data file.  The default is the root name of the input
        fits file appended with an underscore, followed by the
        extension number (ext), followed by the extension ``.txt``.

    cdfile : file path, file object or file-like object, optional
        Output column definitions file.  The default is `None`,
        no column definitions output is produced.

    hfile : file path, file object or file-like object, optional
        Output header parameters file.  The default is `None`,
        no header parameters output is produced.

    ext : int
        The number of the extension containing the table HDU to be
        dumped.

    overwrite : bool, optional
        If ``True``, overwrite the output file if it exists. Raises an
        ``OSError`` if ``False`` and the output file exists. Default is
        ``False``.

        .. versionchanged:: 1.3
           ``overwrite`` replaces the deprecated ``clobber`` argument.

    Notes
    -----
    The primary use for the `tabledump` function is to allow editing in a
    standard text editor of the table data and parameters.  The
    `tableload` function can be used to reassemble the table from the
    three ASCII files.
    """

    # allow file object to already be opened in any of the valid modes
    # and leave the file in the same state (opened or closed) as when
    # the function was called

    mode, closed = _get_file_mode(filename, default='readonly')
    f = fitsopen(filename, mode=mode)

    # Create the default data file name if one was not provided
    try:
        if not datafile:
            root, tail = os.path.splitext(f._file.name)
            datafile = root + '_' + repr(ext) + '.txt'

        # Dump the data from the HDU to the files
        f[ext].dump(datafile, cdfile, hfile, overwrite)
    finally:
        if closed:
            f.close()


if isinstance(tabledump.__doc__, str):
    tabledump.__doc__ += BinTableHDU._tdump_file_format.replace('\n', '\n    ')


def tableload(datafile, cdfile, hfile=None):
    """
    Create a table from the input ASCII files.  The input is from up
    to three separate files, one containing column definitions, one
    containing header parameters, and one containing column data.  The
    header parameters file is not required.  When the header
    parameters file is absent a minimal header is constructed.

    Parameters
    ----------
    datafile : file path, file object or file-like object
        Input data file containing the table data in ASCII format.

    cdfile : file path, file object or file-like object
        Input column definition file containing the names, formats,
        display formats, physical units, multidimensional array
        dimensions, undefined values, scale factors, and offsets
        associated with the columns in the table.

    hfile : file path, file object or file-like object, optional
        Input parameter definition file containing the header
        parameter definitions to be associated with the table.
        If `None`, a minimal header is constructed.

    Notes
    -----
    The primary use for the `tableload` function is to allow the input of
    ASCII data that was edited in a standard text editor of the table
    data and parameters.  The tabledump function can be used to create the
    initial ASCII files.
    """

    return BinTableHDU.load(datafile, cdfile, hfile, replace=True)


if isinstance(tableload.__doc__, str):
    tableload.__doc__ += BinTableHDU._tdump_file_format.replace('\n', '\n    ')


def _getext(filename, mode, *args, ext=None, extname=None, extver=None,
            **kwargs):
    """
    Open the input file, return the `HDUList` and the extension.

    This supports several different styles of extension selection.  See the
    :func:`getdata()` documentation for the different possibilities.
    """

    err_msg = ('Redundant/conflicting extension arguments(s): {}'.format(
            {'args': args, 'ext': ext, 'extname': extname,
             'extver': extver}))

    # This code would be much simpler if just one way of specifying an
    # extension were picked.  But now we need to support all possible ways for
    # the time being.
    if len(args) == 1:
        # Must be either an extension number, an extension name, or an
        # (extname, extver) tuple
        if _is_int(args[0]) or (isinstance(ext, tuple) and len(ext) == 2):
            if ext is not None or extname is not None or extver is not None:
                raise TypeError(err_msg)
            ext = args[0]
        elif isinstance(args[0], str):
            # The first arg is an extension name; it could still be valid
            # to provide an extver kwarg
            if ext is not None or extname is not None:
                raise TypeError(err_msg)
            extname = args[0]
        else:
            # Take whatever we have as the ext argument; we'll validate it
            # below
            ext = args[0]
    elif len(args) == 2:
        # Must be an extname and extver
        if ext is not None or extname is not None or extver is not None:
            raise TypeError(err_msg)
        extname = args[0]
        extver = args[1]
    elif len(args) > 2:
        raise TypeError('Too many positional arguments.')

    if (ext is not None and
            not (_is_int(ext) or
                 (isinstance(ext, tuple) and len(ext) == 2 and
                  isinstance(ext[0], str) and _is_int(ext[1])))):
        raise ValueError(
            'The ext keyword must be either an extension number '
            '(zero-indexed) or a (extname, extver) tuple.')
    if extname is not None and not isinstance(extname, str):
        raise ValueError('The extname argument must be a string.')
    if extver is not None and not _is_int(extver):
        raise ValueError('The extver argument must be an integer.')

    if ext is None and extname is None and extver is None:
        ext = 0
    elif ext is not None and (extname is not None or extver is not None):
        raise TypeError(err_msg)
    elif extname:
        if extver:
            ext = (extname, extver)
        else:
            ext = (extname, 1)
    elif extver and extname is None:
        raise TypeError('extver alone cannot specify an extension.')

    hdulist = fitsopen(filename, mode=mode, **kwargs)

    return hdulist, ext


def _makehdu(data, header):
    if header is None:
        header = Header()
    hdu = _BaseHDU(data, header)
    if hdu.__class__ in (_BaseHDU, _ValidHDU):
        # The HDU type was unrecognized, possibly due to a
        # nonexistent/incomplete header
        if ((isinstance(data, np.ndarray) and data.dtype.fields is not None) or
                isinstance(data, np.recarray)):
            hdu = BinTableHDU(data, header=header)
        elif isinstance(data, np.ndarray):
            hdu = ImageHDU(data, header=header)
        else:
            raise KeyError('Data must be a numpy array.')
    return hdu


def _stat_filename_or_fileobj(filename):
    closed = fileobj_closed(filename)
    name = fileobj_name(filename) or ''

    try:
        loc = filename.tell()
    except AttributeError:
        loc = 0

    noexist_or_empty = ((name and
                         (not os.path.exists(name) or
                          (os.path.getsize(name) == 0)))
                         or (not name and loc == 0))

    return name, closed, noexist_or_empty


def _get_file_mode(filename, default='readonly'):
    """
    Allow file object to already be opened in any of the valid modes and
    and leave the file in the same state (opened or closed) as when
    the function was called.
    """

    mode = default
    closed = fileobj_closed(filename)

    fmode = fileobj_mode(filename)
    if fmode is not None:
        mode = FILE_MODES.get(fmode)
        if mode is None:
            raise OSError(
                "File mode of the input file object ({!r}) cannot be used to "
                "read/write FITS files.".format(fmode))

    return mode, closed
