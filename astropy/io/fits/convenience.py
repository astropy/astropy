import gzip
import os

import numpy as np

from pyfits.file import PYTHON_MODES, _File
from pyfits.hdu.base import _BaseHDU
from pyfits.hdu.hdulist import fitsopen
from pyfits.hdu.image import PrimaryHDU, ImageHDU
from pyfits.hdu.table import BinTableHDU, _TableBaseHDU
from pyfits.header import Header
from pyfits.util import (_with_extensions, deprecated, fileobj_closed,
                         fileobj_name, isfile)


__all__ = ['getheader', 'getdata', 'getval', 'setval', 'delval', 'writeto',
           'append', 'update', 'info', 'tdump', 'tcreate', 'tabledump',
           'tableload']

"""Convenience functions"""


def getheader(filename, *ext, **extkeys):
    """
    Get the header from an extension of a FITS file.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to get header from.  If an opened file object, its mode
        must be one of the following rb, rb+, or ab+).

    ext
        The rest of the arguments are for extension specification.
        `getdata` for explanations/examples.

    Returns
    -------
    header : `Header` object
    """

    mode, closed = _get_file_mode(filename)
    hdulist, _ext = _getext(filename, mode, *ext, **extkeys)
    hdu = hdulist[_ext]
    header = hdu.header

    hdulist.close(closed=closed)
    return header


def getdata(filename, *ext, **extkeys):
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

            >>> getdata('in.fits')

        By extension number::

            >>> getdata('in.fits', 0)    # the primary header
            >>> getdata('in.fits', 2)    # the second extension
            >>> getdata('in.fits', ext=2) # the second extension

        By name, i.e., ``EXTNAME`` value (if unique)::

            >>> getdata('in.fits', 'sci')
            >>> getdata('in.fits', extname='sci') # equivalent

        Note ``EXTNAME`` values are not case sensitive

        By combination of ``EXTNAME`` and EXTVER`` as separate
        arguments or as a tuple::

            >>> getdata('in.fits', 'sci', 2) # EXTNAME='SCI' & EXTVER=2
            >>> getdata('in.fits', extname='sci', extver=2) # equivalent
            >>> getdata('in.fits', ('sci', 2)) # equivalent

        Ambiguous or conflicting specifications will raise an exception::

            >>> getdata('in.fits', ext=('sci',1), extname='err', extver=2)

    lower, upper : bool, optional
        If `lower` or `upper` are `True`, the field names in the
        returned data object will be converted to lower or upper case,
        respectively.

    view : ndarray view class, optional
        When given, the data will be turned wrapped in the given view
        class, by calling::

           data.view(view)

    Returns
    -------
    array : array, record array or groups data object
        Type depends on the type of the extension being referenced.

        If the optional keyword `header` is set to `True`, this
        function will return a (`data`, `header`) tuple.
    """

    if 'header' in extkeys:
        gethdr = extkeys['header']
        del extkeys['header']
    else:
        gethdr = False

    # Code further down rejects unknown keys
    lower = False
    if 'lower' in extkeys:
        lower = extkeys['lower']
        del extkeys['lower']
    upper = False
    if 'upper' in extkeys:
        upper = extkeys['upper']
        del extkeys['upper']
    view = None
    if 'view' in extkeys:
        view = extkeys['view']
        del extkeys['view']

    mode, closed = _get_file_mode(filename)
    hdulist, _ext = _getext(filename, mode, *ext, **extkeys)
    hdu = hdulist[_ext]
    data = hdu.data
    if data is None and _ext == 0:
        try:
            hdu = hdulist[1]
            data = hdu.data
        except IndexError:
            raise IndexError('No data in this HDU.')
    if data is None:
        raise IndexError('No data in this HDU.')
    if gethdr:
        hdr = hdu.header
    hdulist.close(closed=closed)

    # Change case of names if requested
    trans = None
    if lower:
        trans = lambda s: s.lower()
    elif upper:
        trans = lambda s: s.upper()
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
    if view is not None:
        data = data.view(view)

    if gethdr:
        return data, hdr
    else:
        return data


@_with_extensions
def getval(filename, key, *ext, **extkeys):
    """
    Get a keyword's value from a header in a FITS file.

    Parameters
    ----------
    filename : file path, file object, or file like object
        Name of the FITS file, or file object (if opened, mode must be
        one of the following rb, rb+, or ab+).

    key : str
        keyword name

    classExtensions : (optional) **(Deprecated)**
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    ext
        The rest of the arguments are for extension specification.
        See `getdata` for explanations/examples.

    Returns
    -------
    keyword value : string, integer, or float
    """

    hdr = getheader(filename, *ext, **extkeys)
    return hdr[key]


def setval(filename, key, value="", comment=None, before=None, after=None,
           savecomment=False, *ext, **extkeys):
    """
    Set a keyword's value from a header in a FITS file.

    If the keyword already exists, it's value/comment will be updated.
    If it does not exist, a new card will be created and it will be
    placed before or after the specified location.  If no `before` or
    `after` is specified, it will be appended at the end.

    When updating more than one keyword in a file, this convenience
    function is a much less efficient approach compared with opening
    the file for update, modifying the header, and closing the file.

    Parameters
    ----------
    filename : file path, file object, or file like object
        Name of the FITS file, or file object If opened, mode must be
        update (rb+).  An opened file object or `GzipFile` object will
        be closed upon return.

    key : str
        keyword name

    value : str, int, float
        Keyword value, default = ""

    comment : str
        Keyword comment, default = None

    before : str, int
        name of the keyword, or index of the `Card` before which
        the new card will be placed.  The argument `before` takes
        precedence over `after` if both specified. default=`None`.

    after : str, int
        name of the keyword, or index of the `Card` after which the
        new card will be placed. default=`None`.

    savecomment : bool
        when `True`, preserve the current comment for an existing
        keyword.  The argument `savecomment` takes precedence over
        `comment` if both specified.  If `comment` is not specified
        then the current comment will automatically be preserved.
        default=`False`

    ext
        The rest of the arguments are for extension specification.
        See `getdata` for explanations/examples.
    """

    if 'do_not_scale_image_data' not in extkeys:
        extkeys['do_not_scale_image_data'] = True

    hdulist, ext = _getext(filename, mode='update', *ext, **extkeys)
    hdulist[ext].header.update(key, value, comment, before, after, savecomment)
    hdulist.close()


@_with_extensions
def delval(filename, key, *ext, **extkeys):
    """
    Delete all instances of keyword from a header in a FITS file.

    Parameters
    ----------

    filename : file path, file object, or file like object
        Name of the FITS file, or file object If opened, mode must be
        update (rb+).  An opened file object or `GzipFile` object will
        be closed upon return.

    key : str, int
        Keyword name or index

    classExtensions : optional **(Deprecated)**
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    ext
        The rest of the arguments are for extension specification.
        See `getdata` for explanations/examples.
    """

    if 'do_not_scale_image_data' not in extkeys:
        extkeys['do_not_scale_image_data'] = True

    hdulist, ext = _getext(filename, mode='update', *ext, **extkeys)
    del hdulist[ext].header[key]
    hdulist.close()


@_with_extensions
def writeto(filename, data, header=None, **keys):
    """
    Create a new FITS file using the supplied data/header.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened for append (ab+).

    data : array, record array, or groups data object
        data to write to the new file

    header : Header object, optional
        the header associated with `data`. If `None`, a header
        of the appropriate type is created for the supplied data. This
        argument is optional.

    classExtensions : dict, optional **(Deprecated)**
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    clobber : bool, optional
        If `True`, and if filename already exists, it will overwrite
        the file.  Default is `False`.

    checksum : bool, optional
        If `True`, adds both ``DATASUM`` and ``CHECKSUM`` cards to the
        headers of all HDU's written to the file.
    """

    if header is None:
        if 'header' in keys:
            header = keys['header']

    clobber = keys.get('clobber', False)
    output_verify = keys.get('output_verify', 'exception')

    hdu = _makehdu(data, header)
    if not isinstance(hdu, PrimaryHDU) and not isinstance(hdu, _TableBaseHDU):
        hdu = PrimaryHDU(data, header=header)
    checksum = keys.get('checksum', False)
    hdu.writeto(filename, clobber=clobber, output_verify=output_verify,
                checksum=checksum)


@_with_extensions
def append(filename, data, header=None, classExtensions={}, checksum=False,
           verify=True, **keys):
    """
    Append the header/data to FITS file if filename exists, create if not.

    If only `data` is supplied, a minimal header is created.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened for update (rb+)
        unless it is a new file, then it must be opened for append
        (ab+).  A file or `GzipFile` object opened for update will be
        closed after return.

    data : array, table, or group data object
        the new data used for appending

    header : Header object, optional
        The header associated with `data`.  If `None`, an appropriate
        header will be created for the data object supplied.

    classExtensions : dictionary, optional **(Deprecated)**
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    checksum : bool, optional
        When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards to
        the header of the HDU when written to the file.

    verify: bool, optional (True)
        When `True`, the existing FITS file will be read in to verify
        it for correctness before appending.  When `False`, content is
        simply appended to the end of the file.  Setting *verify* to
        `False` can be much faster.
    """

    name, closed, noexist_or_empty = _stat_filename_or_fileobj(filename)

    if noexist_or_empty:
        #
        # The input file or file like object either doesn't exits or is
        # empty.  Use the writeto convenience function to write the
        # output to the empty object.
        #
        writeto(filename, data, header, checksum=checksum, **keys)
    else:
        hdu = _makehdu(data, header)

        if isinstance(hdu, PrimaryHDU):
            hdu = ImageHDU(data, header)

        if verify or not closed:
            f = fitsopen(filename, mode='append')
            f.append(hdu)

            # Set a flag in the HDU so that only this HDU gets a checksum
            # when writing the file.
            hdu._output_checksum = checksum
            f.close(closed=closed)
        else:
            f = _File(filename, mode='append')
            hdu._output_checksum = checksum
            # TODO: Fix this once an API for writing an HDU to a file is
            # settled on
            hdu._writeto(f)
            f.close()


@_with_extensions
def update(filename, data, *ext, **extkeys):
    """
    Update the specified extension with the input data/header.

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to update.  If opened, mode must be update (rb+).  An
        opened file object or `GzipFile` object will be closed upon
        return.

    data : array, table, or group data object
        the new data used for updating

    classExtensions : dict, optional **(Deprecated)**
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    ext
        The rest of the arguments are flexible: the 3rd argument can
        be the header associated with the data.  If the 3rd argument
        is not a `Header`, it (and other positional arguments) are
        assumed to be the extension specification(s).  Header and
        extension specs can also be keyword arguments.  For example::

            >>> update(file, dat, hdr, 'sci')  # update the 'sci' extension
            >>> update(file, dat, 3)  # update the 3rd extension
            >>> update(file, dat, hdr, 3)  # update the 3rd extension
            >>> update(file, dat, 'sci', 2)  # update the 2nd SCI extension
            >>> update(file, dat, 3, header=hdr)  # update the 3rd extension
            >>> update(file, dat, header=hdr, ext=5)  # update the 5th extension
    """

    # parse the arguments
    header = None
    if len(ext) > 0:
        if isinstance(ext[0], Header):
            header = ext[0]
            ext = ext[1:]
        elif not isinstance(ext[0], (int, long, np.integer, str, tuple)):
            raise KeyError('Input argument has wrong data type.')

    if 'header' in extkeys:
        header = extkeys['header']
        del extkeys['header']

    new_hdu = _makehdu(data, header)

    closed = fileobj_closed(filename)

    hdulist, _ext = _getext(filename, 'update', *ext, **extkeys)
    hdulist[_ext] = new_hdu

    hdulist.close(closed=closed)


@_with_extensions
def info(filename, classExtensions={}, output=None, **kwargs):
    """
    Print the summary information on a FITS file.

    This includes the name, type, length of header, data shape and type
    for each extension.

    Parameters
    ----------
    filename : file path, file object, or file like object
        FITS file to obtain info from.  If opened, mode must be one of
        the following: rb, rb+, or ab+.

    classExtensions : dict, optional **(Deprecated)**
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    output : file, optional
        File-like object to output the HDU info to.  Outputs to stdout by
        default.

    kwargs : optional keyword arguments

        - **uint** : bool

            Interpret signed integer data where ``BZERO`` is the
            central value and ``BSCALE == 1`` as unsigned integer
            data.  For example, `int16` data with ``BZERO = 32768``
            and ``BSCALE = 1`` would be treated as `uint16` data.

            Note, for backward compatibility, the kwarg **uint16** may
            be used instead.  The kwarg was renamed when support was
            added for integers of any size.

        - **ignore_missing_end** : bool

            Do not issue an exception when opening a file that is
            missing an ``END`` card in the last header.  Default is
            `True`.
    """

    mode, closed = _get_file_mode(filename, default='copyonwrite')
    # Set the default value for the ignore_missing_end parameter
    if not 'ignore_missing_end' in kwargs:
        kwargs['ignore_missing_end'] = True

    f = fitsopen(filename, mode=mode, **kwargs)
    ret = f.info(output=output)

    if closed:
        f.close()

    return ret


@_with_extensions
def tabledump(filename, datafile=None, cdfile=None, hfile=None, ext=1,
              clobber=False, classExtensions={}):
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

    clobber : bool
        Overwrite the output files if they exist.

    classExtensions : dict **(Deprecated)**
        A dictionary that maps pyfits classes to extensions of those
        classes.  When present in the dictionary, the extension class
        will be constructed in place of the pyfits class.

    Notes
    -----
    The primary use for the `tabledump` function is to allow editing in a
    standard text editor of the table data and parameters.  The
    `tcreate` function can be used to reassemble the table from the
    three ASCII files.
    """

    # allow file object to already be opened in any of the valid modes
    # and leave the file in the same state (opened or closed) as when
    # the function was called

    mode, closed = _get_file_mode(filename, default='copyonwrite')
    f = fitsopen(filename, mode=mode)

    # Create the default data file name if one was not provided

    if not datafile:
        root, tail = os.path.splitext(f._HDUList__file.name)
        datafile = root + '_' + repr(ext) + '.txt'

    # Dump the data from the HDU to the files
    f[ext].dump(datafile, cdfile, hfile, clobber)

    if closed:
        f.close()
tabledump.__doc__ += BinTableHDU.tdump_file_format.replace('\n', '\n    ')
tdump = deprecated(name='tdump')(tabledump)


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
tableload.__doc__ += BinTableHDU.tdump_file_format.replace('\n', '\n    ')
tcreate = deprecated(name='tcreate')(tableload)


def _getext(filename, mode, *ext1, **ext2):
    """Open the input file, return the `HDUList` and the extension."""

    hdulist = fitsopen(filename, mode=mode, **ext2)

    # delete these from the variable keyword argument list so the extension
    # will properly validate
    for key in ['classExtensions', 'ignore_missing_end', 'uint16', 'uint']:
        if key in ext2:
            del ext2[key]

    n_ext1 = len(ext1)
    n_ext2 = len(ext2)

    err_msg = 'Redundant/conflicting keyword arguments(s): %s' % ext2

    # parse the extension spec
    if n_ext1 > 2:
        raise ValueError('too many positional arguments')
    elif n_ext1 == 1:
        if n_ext2 == 0:
            ext = ext1[0]
        else:
            if isinstance(ext1[0], (int, np.integer, tuple)):
                raise KeyError(err_msg)
            if isinstance(ext1[0], basestring):
                if n_ext2 == 1 and 'extver' in ext2:
                    ext = ext1[0], ext2['extver']
                raise KeyError(err_msg)
    elif n_ext1 == 2:
        if n_ext2 == 0:
            ext = ext1
        else:
            raise KeyError(err_msg)
    elif n_ext1 == 0:
        if n_ext2 == 0:
            ext = 0
        elif 'ext' in ext2:
            if n_ext2 == 1:
                ext = ext2['ext']
            elif n_ext2 == 2 and 'extver' in ext2:
                ext = ext2['ext'], ext2['extver']
            else:
                raise KeyError(err_msg)
        else:
            if 'extname' in ext2:
                if 'extver' in ext2:
                    ext = ext2['extname'], ext2['extver']
                else:
                    ext = ext2['extname']
            else:
                raise KeyError('Insufficient keyword arguments: %s' % ext2)

    return hdulist, ext


@_with_extensions
def _makehdu(data, header, classExtensions={}):
    if header is None:
        if ((isinstance(data, np.ndarray) and data.dtype.fields is not None)
            or isinstance(data, np.recarray)):
            hdu = BinTableHDU(data)
        elif isinstance(data, np.ndarray):
            hdu = ImageHDU(data)
        else:
            raise KeyError('Data must be numarray or table data.')
    else:
        hdu = _BaseHDU(data, header)
    return hdu


def _stat_filename_or_fileobj(filename):
    closed = fileobj_closed(filename)
    name = fileobj_name(filename) or ''

    try:
        loc = filename.tell()
    except AttributeError:
        loc = 0

    noexist_or_empty = \
        (name and ((not os.path.exists(name)) or (os.path.getsize(name)==0))) \
         or (not name and loc==0)

    return name, closed, noexist_or_empty


# TODO: Replace this with fileobj_mode
def _get_file_mode(filename, default='readonly'):
    """
    Allow file object to already be opened in any of the valid modes and
    and leave the file in the same state (opened or closed) as when
    the function was called.
    """

    mode = default
    closed = True

    if hasattr(filename, 'closed'):
        closed = filename.closed
    elif hasattr(filename, 'fileobj') and filename.fileobj is not None:
        closed = filename.fileobj.closed

    if (isfile(filename) or
        isinstance(filename, gzip.GzipFile) and not closed):
        if isinstance(filename, gzip.GzipFile):
            file_mode = filename.fileobj.mode
        else:
            file_mode = filename.mode

        for key, val in PYTHON_MODES.iteritems():
            if val == file_mode:
                mode = key
                break

    return mode, closed
