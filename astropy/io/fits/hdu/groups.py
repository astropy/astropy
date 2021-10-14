# Licensed under a 3-clause BSD style license - see PYFITS.rst

import sys
import numpy as np

from .base import DTYPE2BITPIX, DELAYED
from .image import PrimaryHDU
from .table import _TableLikeHDU
from astropy.io.fits.column import Column, ColDefs, FITS2NUMPY
from astropy.io.fits.fitsrec import FITS_rec, FITS_record
from astropy.io.fits.util import _is_int, _is_pseudo_integer, _pseudo_zero

from astropy.utils import lazyproperty


class Group(FITS_record):
    """
    One group of the random group data.
    """

    def __init__(self, input, row=0, start=None, end=None, step=None,
                 base=None):
        super().__init__(input, row, start, end, step, base)

    @property
    def parnames(self):
        return self.array.parnames

    @property
    def data(self):
        # The last column in the coldefs is the data portion of the group
        return self.field(self.array._coldefs.names[-1])

    @lazyproperty
    def _unique(self):
        return _par_indices(self.parnames)

    def par(self, parname):
        """
        Get the group parameter value.
        """

        if _is_int(parname):
            result = self.array[self.row][parname]
        else:
            indx = self._unique[parname.upper()]
            if len(indx) == 1:
                result = self.array[self.row][indx[0]]

            # if more than one group parameter have the same name
            else:
                result = self.array[self.row][indx[0]].astype('f8')
                for i in indx[1:]:
                    result += self.array[self.row][i]

        return result

    def setpar(self, parname, value):
        """
        Set the group parameter value.
        """

        # TODO: It would be nice if, instead of requiring a multi-part value to
        # be an array, there were an *option* to automatically split the value
        # into multiple columns if it doesn't already fit in the array data
        # type.

        if _is_int(parname):
            self.array[self.row][parname] = value
        else:
            indx = self._unique[parname.upper()]
            if len(indx) == 1:
                self.array[self.row][indx[0]] = value

            # if more than one group parameter have the same name, the
            # value must be a list (or tuple) containing arrays
            else:
                if isinstance(value, (list, tuple)) and \
                   len(indx) == len(value):
                    for i in range(len(indx)):
                        self.array[self.row][indx[i]] = value[i]
                else:
                    raise ValueError('Parameter value must be a sequence with '
                                     '{} arrays/numbers.'.format(len(indx)))


class GroupData(FITS_rec):
    """
    Random groups data object.

    Allows structured access to FITS Group data in a manner analogous
    to tables.
    """

    _record_type = Group

    def __new__(cls, input=None, bitpix=None, pardata=None, parnames=[],
                bscale=None, bzero=None, parbscales=None, parbzeros=None):
        """
        Parameters
        ----------
        input : array or FITS_rec instance
            input data, either the group data itself (a
            `numpy.ndarray`) or a record array (`FITS_rec`) which will
            contain both group parameter info and the data.  The rest
            of the arguments are used only for the first case.

        bitpix : int
            data type as expressed in FITS ``BITPIX`` value (8, 16, 32,
            64, -32, or -64)

        pardata : sequence of array
            parameter data, as a list of (numeric) arrays.

        parnames : sequence of str
            list of parameter names.

        bscale : int
            ``BSCALE`` of the data

        bzero : int
            ``BZERO`` of the data

        parbscales : sequence of int
            list of bscales for the parameters

        parbzeros : sequence of int
            list of bzeros for the parameters
        """

        if not isinstance(input, FITS_rec):
            if pardata is None:
                npars = 0
            else:
                npars = len(pardata)

            if parbscales is None:
                parbscales = [None] * npars
            if parbzeros is None:
                parbzeros = [None] * npars

            if parnames is None:
                parnames = [f'PAR{idx + 1}' for idx in range(npars)]

            if len(parnames) != npars:
                raise ValueError('The number of parameter data arrays does '
                                 'not match the number of parameters.')

            unique_parnames = _unique_parnames(parnames + ['DATA'])

            if bitpix is None:
                bitpix = DTYPE2BITPIX[input.dtype.name]

            fits_fmt = GroupsHDU._bitpix2tform[bitpix]  # -32 -> 'E'
            format = FITS2NUMPY[fits_fmt]  # 'E' -> 'f4'
            data_fmt = f'{str(input.shape[1:])}{format}'
            formats = ','.join(([format] * npars) + [data_fmt])
            gcount = input.shape[0]

            cols = [Column(name=unique_parnames[idx], format=fits_fmt,
                           bscale=parbscales[idx], bzero=parbzeros[idx])
                    for idx in range(npars)]
            cols.append(Column(name=unique_parnames[-1], format=fits_fmt,
                               bscale=bscale, bzero=bzero))

            coldefs = ColDefs(cols)

            self = FITS_rec.__new__(cls,
                                    np.rec.array(None,
                                                 formats=formats,
                                                 names=coldefs.names,
                                                 shape=gcount))

            # By default the data field will just be 'DATA', but it may be
            # uniquified if 'DATA' is already used by one of the group names
            self._data_field = unique_parnames[-1]

            self._coldefs = coldefs
            self.parnames = parnames

            for idx, name in enumerate(unique_parnames[:-1]):
                column = coldefs[idx]
                # Note: _get_scale_factors is used here and in other cases
                # below to determine whether the column has non-default
                # scale/zero factors.
                # TODO: Find a better way to do this than using this interface
                scale, zero = self._get_scale_factors(column)[3:5]
                if scale or zero:
                    self._cache_field(name, pardata[idx])
                else:
                    np.rec.recarray.field(self, idx)[:] = pardata[idx]

            column = coldefs[self._data_field]
            scale, zero = self._get_scale_factors(column)[3:5]
            if scale or zero:
                self._cache_field(self._data_field, input)
            else:
                np.rec.recarray.field(self, npars)[:] = input
        else:
            self = FITS_rec.__new__(cls, input)
            self.parnames = None
        return self

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)
        if isinstance(obj, GroupData):
            self.parnames = obj.parnames
        elif isinstance(obj, FITS_rec):
            self.parnames = obj._coldefs.names

    def __getitem__(self, key):
        out = super().__getitem__(key)
        if isinstance(out, GroupData):
            out.parnames = self.parnames
        return out

    @property
    def data(self):
        """
        The raw group data represented as a multi-dimensional `numpy.ndarray`
        array.
        """

        # The last column in the coldefs is the data portion of the group
        return self.field(self._coldefs.names[-1])

    @lazyproperty
    def _unique(self):
        return _par_indices(self.parnames)

    def par(self, parname):
        """
        Get the group parameter values.
        """

        if _is_int(parname):
            result = self.field(parname)
        else:
            indx = self._unique[parname.upper()]
            if len(indx) == 1:
                result = self.field(indx[0])

            # if more than one group parameter have the same name
            else:
                result = self.field(indx[0]).astype('f8')
                for i in indx[1:]:
                    result += self.field(i)

        return result


class GroupsHDU(PrimaryHDU, _TableLikeHDU):
    """
    FITS Random Groups HDU class.

    See the :ref:`astropy:random-groups` section in the Astropy documentation
    for more details on working with this type of HDU.
    """

    _bitpix2tform = {8: 'B', 16: 'I', 32: 'J', 64: 'K', -32: 'E', -64: 'D'}
    _data_type = GroupData
    _data_field = 'DATA'
    """
    The name of the table record array field that will contain the group data
    for each group; 'DATA' by default, but may be preceded by any number of
    underscores if 'DATA' is already a parameter name
    """

    def __init__(self, data=None, header=None):
        super().__init__(data=data, header=header)
        if data is not DELAYED:
            self.update_header()

        # Update the axes; GROUPS HDUs should always have at least one axis
        if len(self._axes) <= 0:
            self._axes = [0]
            self._header['NAXIS'] = 1
            self._header.set('NAXIS1', 0, after='NAXIS')

    @classmethod
    def match_header(cls, header):
        keyword = header.cards[0].keyword
        return (keyword == 'SIMPLE' and 'GROUPS' in header and
                header['GROUPS'] is True)

    @lazyproperty
    def data(self):
        """
        The data of a random group FITS file will be like a binary table's
        data.
        """

        if self._axes == [0]:
            return

        data = self._get_tbdata()
        data._coldefs = self.columns
        data.parnames = self.parnames
        del self.columns
        return data

    @lazyproperty
    def parnames(self):
        """The names of the group parameters as described by the header."""

        pcount = self._header['PCOUNT']
        # The FITS standard doesn't really say what to do if a parname is
        # missing, so for now just assume that won't happen
        return [self._header['PTYPE' + str(idx + 1)] for idx in range(pcount)]

    @lazyproperty
    def columns(self):
        if self._has_data and hasattr(self.data, '_coldefs'):
            return self.data._coldefs

        format = self._bitpix2tform[self._header['BITPIX']]
        pcount = self._header['PCOUNT']
        parnames = []
        bscales = []
        bzeros = []

        for idx in range(pcount):
            bscales.append(self._header.get('PSCAL' + str(idx + 1), None))
            bzeros.append(self._header.get('PZERO' + str(idx + 1), None))
            parnames.append(self._header['PTYPE' + str(idx + 1)])

        formats = [format] * len(parnames)
        dim = [None] * len(parnames)

        # Now create columns from collected parameters, but first add the DATA
        # column too, to contain the group data.
        parnames.append('DATA')
        bscales.append(self._header.get('BSCALE'))
        bzeros.append(self._header.get('BZEROS'))
        data_shape = self.shape[:-1]
        formats.append(str(int(np.prod(data_shape))) + format)
        dim.append(data_shape)
        parnames = _unique_parnames(parnames)

        self._data_field = parnames[-1]

        cols = [Column(name=name, format=fmt, bscale=bscale, bzero=bzero,
                       dim=dim)
                for name, fmt, bscale, bzero, dim in
                zip(parnames, formats, bscales, bzeros, dim)]

        coldefs = ColDefs(cols)
        return coldefs

    @property
    def _nrows(self):
        if not self._data_loaded:
            # The number of 'groups' equates to the number of rows in the table
            # representation of the data
            return self._header.get('GCOUNT', 0)
        else:
            return len(self.data)

    @lazyproperty
    def _theap(self):
        # Only really a lazyproperty for symmetry with _TableBaseHDU
        return 0

    @property
    def is_image(self):
        return False

    @property
    def size(self):
        """
        Returns the size (in bytes) of the HDU's data part.
        """

        size = 0
        naxis = self._header.get('NAXIS', 0)

        # for random group image, NAXIS1 should be 0, so we skip NAXIS1.
        if naxis > 1:
            size = 1
            for idx in range(1, naxis):
                size = size * self._header['NAXIS' + str(idx + 1)]
            bitpix = self._header['BITPIX']
            gcount = self._header.get('GCOUNT', 1)
            pcount = self._header.get('PCOUNT', 0)
            size = abs(bitpix) * gcount * (pcount + size) // 8
        return size

    def update_header(self):
        old_naxis = self._header.get('NAXIS', 0)

        if self._data_loaded:
            if isinstance(self.data, GroupData):
                self._axes = list(self.data.data.shape)[1:]
                self._axes.reverse()
                self._axes = [0] + self._axes
                field0 = self.data.dtype.names[0]
                field0_code = self.data.dtype.fields[field0][0].name
            elif self.data is None:
                self._axes = [0]
                field0_code = 'uint8'  # For lack of a better default
            else:
                raise ValueError('incorrect array type')

            self._header['BITPIX'] = DTYPE2BITPIX[field0_code]

        self._header['NAXIS'] = len(self._axes)

        # add NAXISi if it does not exist
        for idx, axis in enumerate(self._axes):
            if (idx == 0):
                after = 'NAXIS'
            else:
                after = 'NAXIS' + str(idx)

            self._header.set('NAXIS' + str(idx + 1), axis, after=after)

        # delete extra NAXISi's
        for idx in range(len(self._axes) + 1, old_naxis + 1):
            try:
                del self._header['NAXIS' + str(idx)]
            except KeyError:
                pass

        if self._has_data and isinstance(self.data, GroupData):
            self._header.set('GROUPS', True,
                             after='NAXIS' + str(len(self._axes)))
            self._header.set('PCOUNT', len(self.data.parnames), after='GROUPS')
            self._header.set('GCOUNT', len(self.data), after='PCOUNT')

            column = self.data._coldefs[self._data_field]
            scale, zero = self.data._get_scale_factors(column)[3:5]
            if scale:
                self._header.set('BSCALE', column.bscale)
            if zero:
                self._header.set('BZERO', column.bzero)

            for idx, name in enumerate(self.data.parnames):
                self._header.set('PTYPE' + str(idx + 1), name)
                column = self.data._coldefs[idx]
                scale, zero = self.data._get_scale_factors(column)[3:5]
                if scale:
                    self._header.set('PSCAL' + str(idx + 1), column.bscale)
                if zero:
                    self._header.set('PZERO' + str(idx + 1), column.bzero)

        # Update the position of the EXTEND keyword if it already exists
        if 'EXTEND' in self._header:
            if len(self._axes):
                after = 'NAXIS' + str(len(self._axes))
            else:
                after = 'NAXIS'
            self._header.set('EXTEND', after=after)

    def _writedata_internal(self, fileobj):
        """
        Basically copy/pasted from `_ImageBaseHDU._writedata_internal()`, but
        we have to get the data's byte order a different way...

        TODO: Might be nice to store some indication of the data's byte order
        as an attribute or function so that we don't have to do this.
        """

        size = 0

        if self.data is not None:
            self.data._scale_back()

            # Based on the system type, determine the byteorders that
            # would need to be swapped to get to big-endian output
            if sys.byteorder == 'little':
                swap_types = ('<', '=')
            else:
                swap_types = ('<',)
            # deal with unsigned integer 16, 32 and 64 data
            if _is_pseudo_integer(self.data.dtype):
                # Convert the unsigned array to signed
                output = np.array(
                    self.data - _pseudo_zero(self.data.dtype),
                    dtype=f'>i{self.data.dtype.itemsize}')
                should_swap = False
            else:
                output = self.data
                fname = self.data.dtype.names[0]
                byteorder = self.data.dtype.fields[fname][0].str[0]
                should_swap = (byteorder in swap_types)

            if should_swap:
                if output.flags.writeable:
                    output.byteswap(True)
                    try:
                        fileobj.writearray(output)
                    finally:
                        output.byteswap(True)
                else:
                    # For read-only arrays, there is no way around making
                    # a byteswapped copy of the data.
                    fileobj.writearray(output.byteswap(False))
            else:
                fileobj.writearray(output)

            size += output.size * output.itemsize
        return size

    def _verify(self, option='warn'):
        errs = super()._verify(option=option)

        # Verify locations and values of mandatory keywords.
        self.req_cards('NAXIS', 2,
                       lambda v: (_is_int(v) and 1 <= v <= 999), 1,
                       option, errs)
        self.req_cards('NAXIS1', 3, lambda v: (_is_int(v) and v == 0), 0,
                       option, errs)

        after = self._header['NAXIS'] + 3
        pos = lambda x: x >= after

        self.req_cards('GCOUNT', pos, _is_int, 1, option, errs)
        self.req_cards('PCOUNT', pos, _is_int, 0, option, errs)
        self.req_cards('GROUPS', pos, lambda v: (v is True), True, option,
                       errs)
        return errs

    def _calculate_datasum(self):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._has_data:

            # We have the data to be used.

            # Check the byte order of the data.  If it is little endian we
            # must swap it before calculating the datasum.
            # TODO: Maybe check this on a per-field basis instead of assuming
            # that all fields have the same byte order?
            byteorder = \
                self.data.dtype.fields[self.data.dtype.names[0]][0].str[0]

            if byteorder != '>':
                if self.data.flags.writeable:
                    byteswapped = True
                    d = self.data.byteswap(True)
                    d.dtype = d.dtype.newbyteorder('>')
                else:
                    # If the data is not writeable, we just make a byteswapped
                    # copy and don't bother changing it back after
                    d = self.data.byteswap(False)
                    d.dtype = d.dtype.newbyteorder('>')
                    byteswapped = False
            else:
                byteswapped = False
                d = self.data

            byte_data = d.view(type=np.ndarray, dtype=np.ubyte)

            cs = self._compute_checksum(byte_data)

            # If the data was byteswapped in this method then return it to
            # its original little-endian order.
            if byteswapped:
                d.byteswap(True)
                d.dtype = d.dtype.newbyteorder('<')

            return cs
        else:
            # This is the case where the data has not been read from the file
            # yet.  We can handle that in a generic manner so we do it in the
            # base class.  The other possibility is that there is no data at
            # all.  This can also be handled in a generic manner.
            return super()._calculate_datasum()

    def _summary(self):
        summary = super()._summary()
        name, ver, classname, length, shape, format, gcount = summary

        # Drop the first axis from the shape
        if shape:
            shape = shape[1:]

            if shape and all(shape):
                # Update the format
                format = self.columns[0].dtype.name

        # Update the GCOUNT report
        gcount = f'{self._gcount} Groups  {self._pcount} Parameters'
        return (name, ver, classname, length, shape, format, gcount)


def _par_indices(names):
    """
    Given a list of objects, returns a mapping of objects in that list to the
    index or indices at which that object was found in the list.
    """

    unique = {}
    for idx, name in enumerate(names):
        # Case insensitive
        name = name.upper()
        if name in unique:
            unique[name].append(idx)
        else:
            unique[name] = [idx]
    return unique


def _unique_parnames(names):
    """
    Given a list of parnames, including possible duplicates, returns a new list
    of parnames with duplicates prepended by one or more underscores to make
    them unique.  This is also case insensitive.
    """

    upper_names = set()
    unique_names = []

    for name in names:
        name_upper = name.upper()
        while name_upper in upper_names:
            name = '_' + name
            name_upper = '_' + name_upper

        unique_names.append(name)
        upper_names.add(name_upper)

    return unique_names
