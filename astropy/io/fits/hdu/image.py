# Licensed under a 3-clause BSD style license - see PYFITS.rst

import sys
import warnings

import numpy as np

from .base import DELAYED, _ValidHDU, ExtensionHDU
from ..header import Header
from ..util import _is_pseudo_unsigned, _unsigned_zero, _is_int
from ..verify import VerifyWarning

from ....extern.six import string_types
from ....utils import isiterable, lazyproperty


class _ImageBaseHDU(_ValidHDU):
    """FITS image HDU base class.

    Attributes
    ----------
    header
        image header

    data
        image data
    """

    # mappings between FITS and numpy typecodes
    # TODO: Maybe make these module-level constants instead...
    NumCode = {8: 'uint8', 16: 'int16', 32: 'int32', 64: 'int64',
               -32: 'float32', -64: 'float64'}
    ImgCode = {'uint8': 8, 'int16': 16, 'uint16': 16, 'int32': 32,
               'uint32': 32, 'int64': 64, 'uint64': 64, 'float32': -32,
               'float64': -64}

    standard_keyword_comments = {
        'SIMPLE': 'conforms to FITS standard',
        'XTENSION': 'Image extension',
        'BITPIX': 'array data type',
        'NAXIS': 'number of array dimensions',
        'GROUPS': 'has groups',
        'PCOUNT': 'number of parameters',
        'GCOUNT': 'number of groups'
    }

    def __init__(self, data=None, header=None, do_not_scale_image_data=False,
                 uint=False, scale_back=False, ignore_blank=False, **kwargs):

        from .groups import GroupsHDU

        super(_ImageBaseHDU, self).__init__(data=data, header=header)

        if header is not None:
            if not isinstance(header, Header):
                # TODO: Instead maybe try initializing a new Header object from
                # whatever is passed in as the header--there are various types
                # of objects that could work for this...
                raise ValueError('header must be a Header object')

        if data is DELAYED:
            # Presumably if data is DELAYED then this HDU is coming from an
            # open file, and was not created in memory
            if header is None:
                # this should never happen
                raise ValueError('No header to setup HDU.')

            # if the file is read the first time, no need to copy, and keep it
            # unchanged
            else:
                self._header = header
        else:
            # TODO: Some of this card manipulation should go into the
            # PrimaryHDU and GroupsHDU subclasses
            # construct a list of cards of minimal header
            if isinstance(self, ExtensionHDU):
                c0 = ('XTENSION', 'IMAGE',
                      self.standard_keyword_comments['XTENSION'])
            else:
                c0 = ('SIMPLE', True, self.standard_keyword_comments['SIMPLE'])
            cards = [
                c0,
                ('BITPIX',    8, self.standard_keyword_comments['BITPIX']),
                ('NAXIS',     0, self.standard_keyword_comments['NAXIS'])]

            if isinstance(self, GroupsHDU):
                cards.append(('GROUPS', True,
                             self.standard_keyword_comments['GROUPS']))

            if isinstance(self, (ExtensionHDU, GroupsHDU)):
                cards.append(('PCOUNT',    0,
                              self.standard_keyword_comments['PCOUNT']))
                cards.append(('GCOUNT',    1,
                              self.standard_keyword_comments['GCOUNT']))

            if header is not None:
                orig = header.copy()
                header = Header(cards)
                header.extend(orig, strip=True, update=True, end=True)
            else:
                header = Header(cards)

            self._header = header

        self._do_not_scale_image_data = do_not_scale_image_data

        self._uint = uint
        self._scale_back = scale_back

        if do_not_scale_image_data:
            self._bzero = 0
            self._bscale = 1
        else:
            self._bzero = self._header.get('BZERO', 0)
            self._bscale = self._header.get('BSCALE', 1)

        # Save off other important values from the header needed to interpret
        # the image data
        self._axes = [self._header.get('NAXIS' + str(axis + 1), 0)
                      for axis in range(self._header.get('NAXIS', 0))]
        self._bitpix = self._header.get('BITPIX', 8)
        self._gcount = self._header.get('GCOUNT', 1)
        self._pcount = self._header.get('PCOUNT', 0)
        self._blank = None if ignore_blank else self._header.get('BLANK')
        self._verify_blank()

        self._orig_bitpix = self._bitpix
        self._orig_bzero = self._bzero
        self._orig_bscale = self._bscale
        self._orig_blank = self._header.get('BLANK')

        # Set the name attribute if it was provided (if this is an ImageHDU
        # this will result in setting the EXTNAME keyword of the header as
        # well)
        if 'name' in kwargs and kwargs['name']:
            self.name = kwargs['name']

        # Set to True if the data or header is replaced, indicating that
        # update_header should be called
        self._modified = False

        if data is DELAYED:
            if (not do_not_scale_image_data and
                    (self._bscale != 1 or self._bzero != 0)):
                # This indicates that when the data is accessed or written out
                # to a new file it will need to be rescaled
                self._data_needs_rescale = True
            return
        else:
            self.data = data
            self.update_header()

    @classmethod
    def match_header(cls, header):
        """
        _ImageBaseHDU is sort of an abstract class for HDUs containing image
        data (as opposed to table data) and should never be used directly.
        """

        raise NotImplementedError

    @property
    def is_image(self):
        return True

    @property
    def section(self):
        """
        Access a section of the image array without loading the entire array
        into memory.  The :class:`Section` object returned by this attribute is
        not meant to be used directly by itself.  Rather, slices of the section
        return the appropriate slice of the data, and loads *only* that section
        into memory.

        Sections are mostly obsoleted by memmap support, but should still be
        used to deal with very large scaled images.  See the
        :ref:`data-sections` section of the PyFITS documentation for more
        details.
        """

        return Section(self)

    @property
    def shape(self):
        """
        Shape of the image array--should be equivalent to ``self.data.shape``.
        """

        # Determine from the values read from the header
        return tuple(reversed(self._axes))

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, header):
        self._header = header
        self._modified = True
        self.update_header()

    @lazyproperty
    def data(self):
        """
        Image/array data as a `~numpy.ndarray`.

        Please remember that the order of axes on an Numpy array are opposite
        of the order specified in the FITS file.  For example for a 2D image
        the "rows" or y-axis are the first dimension, and the "columns" or
        x-axis are the second dimension.

        If the data is scaled using the BZERO and BSCALE parameters, this
        attribute returns the data scaled to its physical values unless the
        file was opened with ``do_not_scale_image_data=True``.
        """

        if len(self._axes) < 1:
            return

        data = self._get_scaled_image_data(self._data_offset, self.shape)
        self._update_header_scale_info(data.dtype)

        return data

    @data.setter
    def data(self, data):
        if 'data' in self.__dict__:
            if self.__dict__['data'] is data:
                return
            else:
                self._data_replaced = True
        else:
            self._data_replaced = True

        if data is not None and not isinstance(data, np.ndarray):
            # Try to coerce the data into a numpy array--this will work, on
            # some level, for most objects
            try:
                data = np.array(data)
            except:
                raise TypeError('data object %r could not be coerced into an '
                                'ndarray' % data)

        self.__dict__['data'] = data
        self._modified = True

        if isinstance(data, np.ndarray):
            self._bitpix = _ImageBaseHDU.ImgCode[data.dtype.name]
            self._orig_bitpix = self._bitpix
            self._orig_bscale = 1
            self._orig_bzero = 0
            self._axes = list(data.shape)
            self._axes.reverse()
        elif self.data is None:
            self._axes = []
        else:
            raise ValueError('not a valid data array')

        self.update_header()

        # returning the data signals to lazyproperty that we've already handled
        # setting self.__dict__['data']
        return data

    def update_header(self):
        """
        Update the header keywords to agree with the data.
        """

        if not (self._modified or self._header._modified or
                (self._has_data and self.shape != self.data.shape)):
            # Not likely that anything needs updating
            return

        old_naxis = self._header.get('NAXIS', 0)

        if 'BITPIX' not in self._header:
            bitpix_comment = self.standard_keyword_comments['BITPIX']
        else:
            bitpix_comment = self._header.comments['BITPIX']

        # Update the BITPIX keyword and ensure it's in the correct
        # location in the header
        self._header.set('BITPIX', self._bitpix, bitpix_comment, after=0)

        # If the data's shape has changed (this may have happened without our
        # noticing either via a direct update to the data.shape attribute) we
        # need to update the internal self._axes
        if self._has_data and self.shape != self.data.shape:
            self._axes = list(self.data.shape)
            self._axes.reverse()

        # Update the NAXIS keyword and ensure it's in the correct location in
        # the header
        if 'NAXIS' in self._header:
            naxis_comment = self._header.comments['NAXIS']
        else:
            naxis_comment = self.standard_keyword_comments['NAXIS']
        self._header.set('NAXIS', len(self._axes), naxis_comment,
                         after='BITPIX')

        # TODO: This routine is repeated in several different classes--it
        # should probably be made available as a method on all standard HDU
        # types
        # add NAXISi if it does not exist
        for idx, axis in enumerate(self._axes):
            naxisn = 'NAXIS' + str(idx + 1)
            if naxisn in self._header:
                self._header[naxisn] = axis
            else:
                if (idx == 0):
                    after = 'NAXIS'
                else:
                    after = 'NAXIS' + str(idx)
                self._header.set(naxisn, axis, after=after)

        # delete extra NAXISi's
        for idx in range(len(self._axes) + 1, old_naxis + 1):
            try:
                del self._header['NAXIS' + str(idx)]
            except KeyError:
                pass

        if 'BLANK' in self._header:
            self._blank = self._header['BLANK']

        self._update_uint_scale_keywords()

        self._modified = False

    def _update_header_scale_info(self, dtype=None):
        if (not self._do_not_scale_image_data and
                not (self._orig_bzero == 0 and self._orig_bscale == 1)):

            if dtype is None:
                dtype = self._dtype_for_bitpix()

            if (dtype is not None and dtype.kind == 'u' and
                    (self._scale_back or self._scale_back is None)):
                # Data is pseudo-unsigned integers, and the scale_back option
                # was not explicitly set to False, so preserve all the scale
                # factors
                return

            for keyword in ['BSCALE', 'BZERO']:
                try:
                    del self._header[keyword]
                    # Since _update_header_scale_info can, currently, be called
                    # *after* _prewriteto(), replace these with blank cards so
                    # the header size doesn't change
                    self._header.append()
                except KeyError:
                    pass

            if dtype is None:
                dtype = self._dtype_for_bitpix()
            if dtype is not None:
                self._header['BITPIX'] = _ImageBaseHDU.ImgCode[dtype.name]

            self._bzero = 0
            self._bscale = 1
            self._bitpix = self._header['BITPIX']
            self._blank = self._header.pop('BLANK', None)

    def scale(self, type=None, option='old', bscale=1, bzero=0):
        """
        Scale image data by using ``BSCALE``/``BZERO``.

        Call to this method will scale `data` and update the keywords of
        ``BSCALE`` and ``BZERO`` in the HDU's header.  This method should only
        be used right before writing to the output file, as the data will be
        scaled and is therefore not very usable after the call.

        Parameters
        ----------
        type : str, optional
            destination data type, use a string representing a numpy
            dtype name, (e.g. ``'uint8'``, ``'int16'``, ``'float32'``
            etc.).  If is `None`, use the current data type.

        option : str
            How to scale the data: if ``"old"``, use the original
            ``BSCALE`` and ``BZERO`` values when the data was
            read/created. If ``"minmax"``, use the minimum and maximum
            of the data to scale.  The option will be overwritten by
            any user specified ``bscale``/``bzero`` values.

        bscale, bzero : int, optional
            User-specified ``BSCALE`` and ``BZERO`` values
        """

        # Disable blank support for now
        self._scale_internal(type=type, option=option, bscale=bscale,
                             bzero=bzero, blank=None)

    def _scale_internal(self, type=None, option='old', bscale=1, bzero=0,
                        blank=0):
        """
        This is an internal implementation of the `scale` method, which
        also supports handling BLANK properly.

        TODO: This is only needed for fixing #3865 without introducing any
        public API changes.  We should support BLANK better when rescaling
        data, and when that is added the need for this internal interface
        should go away.

        Note: the default of ``blank=0`` merely reflects the current behavior,
        and is not necessarily a deliberate choice (better would be to disallow
        conversion of floats to ints without specifying a BLANK if there are
        NaN/inf values).
        """

        if self.data is None:
            return

        # Determine the destination (numpy) data type
        if type is None:
            type = self.NumCode[self._bitpix]
        _type = getattr(np, type)

        # Determine how to scale the data
        # bscale and bzero takes priority
        if (bscale != 1 or bzero != 0):
            _scale = bscale
            _zero = bzero
        else:
            if option == 'old':
                _scale = self._orig_bscale
                _zero = self._orig_bzero
            elif option == 'minmax':
                if issubclass(_type, np.floating):
                    _scale = 1
                    _zero = 0
                else:

                    min = np.minimum.reduce(self.data.flat)
                    max = np.maximum.reduce(self.data.flat)

                    if _type == np.uint8:  # uint8 case
                        _zero = min
                        _scale = (max - min) / (2.0 ** 8 - 1)
                    else:
                        _zero = (max + min) / 2.0

                        # throw away -2^N
                        nbytes = 8 * _type().itemsize
                        _scale = (max - min) / (2.0 ** nbytes - 2)

        # Do the scaling
        if _zero != 0:
            # 0.9.6.3 to avoid out of range error for BZERO = +32768
            self.data += -_zero
            self._header['BZERO'] = _zero
        else:
            try:
                del self._header['BZERO']
            except KeyError:
                pass

        if _scale and _scale != 1:
            self.data = self.data / _scale
            self._header['BSCALE'] = _scale
        else:
            try:
                del self._header['BSCALE']
            except KeyError:
                pass

        # Set blanks
        if blank is not None and issubclass(_type, np.integer):
            # TODO: Perhaps check that the requested BLANK value fits in the
            # integer type being scaled to?
            self.data[np.isnan(self.data)] = blank
            self._header['BLANK'] = blank

        if self.data.dtype.type != _type:
            self.data = np.array(np.around(self.data), dtype=_type)

        # Update the BITPIX Card to match the data
        self._bitpix = _ImageBaseHDU.ImgCode[self.data.dtype.name]
        self._bzero = self._header.get('BZERO', 0)
        self._bscale = self._header.get('BSCALE', 1)
        self._blank = blank
        self._header['BITPIX'] = self._bitpix

        # Since the image has been manually scaled, the current
        # bitpix/bzero/bscale now serve as the 'original' scaling of the image,
        # as though the original image has been completely replaced
        self._orig_bitpix = self._bitpix
        self._orig_bzero = self._bzero
        self._orig_bscale = self._bscale
        self._orig_blank = self._blank

    def _verify(self, option='warn'):
        # update_header can fix some things that would otherwise cause
        # verification to fail, so do that now...
        self.update_header()
        self._verify_blank()

        return super(_ImageBaseHDU, self)._verify(option)

    def _verify_blank(self):
        # Probably not the best place for this (it should probably happen
        # in _verify as well) but I want to be able to raise this warning
        # both when the HDU is created and when written
        if self._blank is None:
            return

        messages = []
        # TODO: Once the FITSSchema framewhere is merged these warnings
        # should be handled by the schema
        if not _is_int(self._blank):
            messages.append(
                "Invalid value for 'BLANK' keyword in header: {0!r} "
                "The 'BLANK' keyword must be an integer.  It will be "
                "ignored in the meantime.".format(self._blank))
            self._blank = None
        if not self._bitpix > 0:
            messages.append(
                "Invalid 'BLANK' keyword in header.  The 'BLANK' keyword "
                "is only applicable to integer data, and will be ignored "
                "in this HDU.")
            self._blank = None

        for msg in messages:
            warnings.warn(msg, VerifyWarning)

    def _prewriteto(self, checksum=False, inplace=False):
        if self._scale_back:
            self._scale_internal(self.NumCode[self._orig_bitpix],
                                 blank=self._orig_blank)

        self.update_header()
        if not inplace and self._data_needs_rescale:
            # Go ahead and load the scaled image data and update the header
            # with the correct post-rescaling headers
            _ = self.data

        return super(_ImageBaseHDU, self)._prewriteto(checksum, inplace)

    def _writedata_internal(self, fileobj):
        size = 0

        if self.data is not None:
            # Based on the system type, determine the byteorders that
            # would need to be swapped to get to big-endian output
            if sys.byteorder == 'little':
                swap_types = ('<', '=')
            else:
                swap_types = ('<',)
            # deal with unsigned integer 16, 32 and 64 data
            if _is_pseudo_unsigned(self.data.dtype):
                # Convert the unsigned array to signed
                output = np.array(
                    self.data - _unsigned_zero(self.data.dtype),
                    dtype='>i%d' % self.data.dtype.itemsize)
                should_swap = False
            else:
                output = self.data
                byteorder = output.dtype.str[0]
                should_swap = (byteorder in swap_types)

            if not fileobj.simulateonly:
                if should_swap:
                    output.byteswap(True)
                    try:
                        fileobj.writearray(output)
                    finally:
                        output.byteswap(True)
                else:
                    fileobj.writearray(output)

            size += output.size * output.itemsize

        return size

    def _dtype_for_bitpix(self):
        """
        Determine the dtype that the data should be converted to depending on
        the BITPIX value in the header, and possibly on the BSCALE value as
        well.  Returns None if there should not be any change.
        """

        bitpix = self._orig_bitpix
        # Handle possible conversion to uints if enabled
        if self._uint and self._orig_bscale == 1:
            for bits, dtype in ((16, np.dtype('uint16')),
                                (32, np.dtype('uint32')),
                                (64, np.dtype('uint64'))):
                if bitpix == bits and self._orig_bzero == 1 << (bits - 1):
                    return dtype

        if bitpix > 16:  # scale integers to Float64
            return np.dtype('float64')
        elif bitpix > 0:  # scale integers to Float32
            return np.dtype('float32')

    def _convert_pseudo_unsigned(self, data):
        """
        Handle "pseudo-unsigned" integers, if the user requested it.  Returns
        the converted data array if so; otherwise returns None.

        In this case case, we don't need to handle BLANK to convert it to NAN,
        since we can't do NaNs with integers, anyway, i.e. the user is
        responsible for managing blanks.
        """

        dtype = self._dtype_for_bitpix()
        # bool(dtype) is always False--have to explicitly compare to None; this
        # caused a fair amount of hair loss
        if dtype is not None and dtype.kind == 'u':
            # Convert the input raw data into an unsigned integer array and
            # then scale the data adjusting for the value of BZERO.  Note that
            # we subtract the value of BZERO instead of adding because of the
            # way numpy converts the raw signed array into an unsigned array.
            bits = dtype.itemsize * 8
            data = np.array(data, dtype=dtype)
            data -= np.uint64(1 << (bits - 1))

            return data

    def _get_scaled_image_data(self, offset, shape):
        """
        Internal function for reading image data from a file and apply scale
        factors to it.  Normally this is used for the entire image, but it
        supports alternate offset/shape for Section support.
        """

        code = _ImageBaseHDU.NumCode[self._orig_bitpix]

        raw_data = self._get_raw_data(shape, code, offset)
        raw_data.dtype = raw_data.dtype.newbyteorder('>')

        if self._do_not_scale_image_data or (
                self._orig_bzero == 0 and self._orig_bscale == 1 and
                self._blank is None):
            # No further conversion of the data is necessary
            return raw_data

        try:
            if self._file.strict_memmap:
                raise ValueError("Cannot load a memory-mapped image: "
                                 "BZERO/BSCALE/BLANK header keywords present. "
                                 "Set memmap=False.")
        except AttributeError:  # strict_memmap not set
            pass

        data = None
        if not (self._orig_bzero == 0 and self._orig_bscale == 1):
            data = self._convert_pseudo_unsigned(raw_data)

        if data is None:
            # In these cases, we end up with floating-point arrays and have to
            # apply bscale and bzero. We may have to handle BLANK and convert
            # to NaN in the resulting floating-point arrays.
            # The BLANK keyword should only be applied for integer data (this
            # is checked in __init__ but it can't hurt to double check here)
            blanks = None

            if self._blank is not None and self._bitpix > 0:
                blanks = raw_data.flat == self._blank
                # The size of blanks in bytes is the number of elements in
                # raw_data.flat.  However, if we use np.where instead we will
                # only use 8 bytes for each index where the condition is true.
                # So if the number of blank items is fewer than
                # len(raw_data.flat) / 8, using np.where will use less memory
                if blanks.sum() < len(blanks) / 8:
                    blanks = np.where(blanks)

            new_dtype = self._dtype_for_bitpix()
            if new_dtype is not None:
                data = np.array(raw_data, dtype=new_dtype)
            else:  # floating point cases
                if self._file is not None and self._file.memmap:
                    data = raw_data.copy()
                elif not raw_data.flags.writeable:
                    # create a writeable copy if needed
                    data = raw_data.copy()
                # if not memmap, use the space already in memory
                else:
                    data = raw_data

            del raw_data

            if self._orig_bscale != 1:
                np.multiply(data, self._orig_bscale, data)
            if self._orig_bzero != 0:
                data += self._orig_bzero

            if self._blank:
                data.flat[blanks] = np.nan

        return data

    # TODO: Move the GroupsHDU-specific summary code to GroupsHDU itself
    def _summary(self):
        """
        Summarize the HDU: name, dimensions, and formats.
        """

        class_name = self.__class__.__name__

        # if data is touched, use data info.
        if self._data_loaded:
            if self.data is None:
                format = ''
            else:
                format = self.data.dtype.name
                format = format[format.rfind('.')+1:]
        else:
            if self.shape and all(self.shape):
                # Only show the format if all the dimensions are non-zero
                # if data is not touched yet, use header info.
                format = self.NumCode[self._bitpix]
            else:
                format = ''

        # Display shape in FITS-order
        shape = tuple(reversed(self.shape))

        return (self.name, class_name, len(self._header), shape, format, '')

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._has_data:
            # We have the data to be used.
            d = self.data

            # First handle the special case where the data is unsigned integer
            # 16, 32 or 64
            if _is_pseudo_unsigned(self.data.dtype):
                d = np.array(self.data - _unsigned_zero(self.data.dtype),
                             dtype='i%d' % self.data.dtype.itemsize)

            # Check the byte order of the data.  If it is little endian we
            # must swap it before calculating the datasum.
            if d.dtype.str[0] != '>':
                byteswapped = True
                d = d.byteswap(True)
                d.dtype = d.dtype.newbyteorder('>')
            else:
                byteswapped = False

            cs = self._compute_checksum(d.flatten().view(np.uint8),
                                        blocking=blocking)

            # If the data was byteswapped in this method then return it to
            # its original little-endian order.
            if byteswapped and not _is_pseudo_unsigned(self.data.dtype):
                d.byteswap(True)
                d.dtype = d.dtype.newbyteorder('<')

            return cs
        else:
            # This is the case where the data has not been read from the file
            # yet.  We can handle that in a generic manner so we do it in the
            # base class.  The other possibility is that there is no data at
            # all.  This can also be handled in a generic manner.
            return super(_ImageBaseHDU, self)._calculate_datasum(
                blocking=blocking)


class Section(object):
    """
    Image section.

    Slices of this object load the corresponding section of an image array from
    the underlying FITS file on disk, and applies any BSCALE/BZERO factors.

    Section slices cannot be assigned to, and modifications to a section are
    not saved back to the underlying file.

    See the :ref:`data-sections` section of the PyFITS documentation for more
    details.
    """

    def __init__(self, hdu):
        self.hdu = hdu

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        naxis = len(self.hdu.shape)
        return_scalar = (all(isinstance(k, (int, np.integer)) for k in key)
                         and len(key) == naxis)
        if not any(k is Ellipsis for k in key):
            # We can always add a ... at the end, after making note of whether
            # to return a scalar.
            key += Ellipsis,
        ellipsis_count = len([k for k in key if k is Ellipsis])
        if len(key) - ellipsis_count > naxis or ellipsis_count > 1:
            raise IndexError('too many indices for array')
        # Insert extra dimensions as needed.
        idx = next(i for i, k in enumerate(key + (Ellipsis,)) if k is Ellipsis)
        key = key[:idx] + (slice(None),) * (naxis - len(key) + 1) + key[idx+1:]
        return_0dim = (all(isinstance(k, (int, np.integer)) for k in key)
                       and len(key) == naxis)

        dims = []
        offset = 0
        # Find all leading axes for which a single point is used.
        for idx in range(naxis):
            axis = self.hdu.shape[idx]
            indx = _IndexInfo(key[idx], axis)
            offset = offset * axis + indx.offset
            if not _is_int(key[idx]):
                dims.append(indx.npts)
                break

        is_contiguous = indx.contiguous
        for jdx in range(idx + 1, naxis):
            axis = self.hdu.shape[jdx]
            indx = _IndexInfo(key[jdx], axis)
            dims.append(indx.npts)
            if indx.npts == axis and indx.contiguous:
                # The offset needs to multiply the length of all remaining axes
                offset *= axis
            else:
                is_contiguous = False

        if is_contiguous:
            dims = tuple(dims) or (1,)
            bitpix = self.hdu._orig_bitpix
            offset = self.hdu._data_offset + offset * abs(bitpix) // 8
            data = self.hdu._get_scaled_image_data(offset, dims)
        else:
            data = self._getdata(key)

        if return_scalar:
            data = data.item()
        elif return_0dim:
            data = data.squeeze()
        return data

    def _getdata(self, keys):
        for idx, (key, axis) in enumerate(zip(keys, self.hdu.shape)):
            if isinstance(key, slice):
                ks = range(*key.indices(axis))
                break
            elif isiterable(key):
                # Handle both integer and boolean arrays.
                ks = np.arange(axis, dtype=int)[key]
                break
            # This should always break at some point if _getdata is called.

        data = [self[keys[:idx] + (k,) + keys[idx + 1:]] for k in ks]

        if any(isinstance(key, slice) or isiterable(key)
               for key in keys[idx + 1:]):
            # data contains multidimensional arrays; combine them.
            return np.array(data)
        else:
            # Only singleton dimensions remain; concatenate in a 1D array.
            return np.concatenate([np.atleast_1d(array) for array in data])


class PrimaryHDU(_ImageBaseHDU):
    """
    FITS primary HDU class.
    """

    _default_name = 'PRIMARY'

    def __init__(self, data=None, header=None, do_not_scale_image_data=False,
                 ignore_blank=False,
                 uint=False, scale_back=None):
        """
        Construct a primary HDU.

        Parameters
        ----------
        data : array or DELAYED, optional
            The data in the HDU.

        header : Header instance, optional
            The header to be used (as a template).  If ``header`` is `None`, a
            minimal header will be provided.

        do_not_scale_image_data : bool, optional
            If `True`, image data is not scaled using BSCALE/BZERO values
            when read.

        ignore_blank : bool, optional
            If `True`, the BLANK header keyword will be ignored if present.
            Otherwise, pixels equal to this value will be replaced with
            NaNs.

        uint : bool, optional
            Interpret signed integer data where ``BZERO`` is the
            central value and ``BSCALE == 1`` as unsigned integer
            data.  For example, ``int16`` data with ``BZERO = 32768``
            and ``BSCALE = 1`` would be treated as ``uint16`` data.

        scale_back : bool, optional
            If `True`, when saving changes to a file that contained scaled
            image data, restore the data to the original type and reapply the
            original BSCALE/BZERO values.  This could lead to loss of accuracy
            if scaling back to integer values after performing floating point
            operations on the data.
        """

        super(PrimaryHDU, self).__init__(
            data=data, header=header,
            do_not_scale_image_data=do_not_scale_image_data, uint=uint,
            ignore_blank=ignore_blank,
            scale_back=scale_back)

        # insert the keywords EXTEND
        if header is None:
            dim = self._header['NAXIS']
            if dim == 0:
                dim = ''
            self._header.set('EXTEND', True, after='NAXIS' + str(dim))

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        return (card.keyword == 'SIMPLE' and
                ('GROUPS' not in header or header['GROUPS'] != True) and
                card.value == True)

    def update_header(self):
        super(PrimaryHDU, self).update_header()

        # Update the position of the EXTEND keyword if it already exists
        if 'EXTEND' in self._header:
            if len(self._axes):
                after = 'NAXIS' + str(len(self._axes))
            else:
                after = 'NAXIS'
            self._header.set('EXTEND', after=after)

    def _verify(self, option='warn'):
        errs = super(PrimaryHDU, self)._verify(option=option)

        # Verify location and value of mandatory keywords.
        # The EXTEND keyword is only mandatory if the HDU has extensions; this
        # condition is checked by the HDUList object.  However, if we already
        # have an EXTEND keyword check that its position is correct
        if 'EXTEND' in self._header:
            naxis = self._header.get('NAXIS', 0)
            self.req_cards('EXTEND', naxis + 3, lambda v: isinstance(v, bool),
                           True, option, errs)
        return errs


class ImageHDU(_ImageBaseHDU, ExtensionHDU):
    """
    FITS image extension HDU class.
    """

    _extension = 'IMAGE'

    def __init__(self, data=None, header=None, name=None,
                 do_not_scale_image_data=False, uint=False, scale_back=None):
        """
        Construct an image HDU.

        Parameters
        ----------
        data : array
            The data in the HDU.

        header : Header instance
            The header to be used (as a template).  If ``header`` is
            `None`, a minimal header will be provided.

        name : str, optional
            The name of the HDU, will be the value of the keyword
            ``EXTNAME``.

        do_not_scale_image_data : bool, optional
            If `True`, image data is not scaled using BSCALE/BZERO values
            when read.

        uint : bool, optional
            Interpret signed integer data where ``BZERO`` is the
            central value and ``BSCALE == 1`` as unsigned integer
            data.  For example, ``int16`` data with ``BZERO = 32768``
            and ``BSCALE = 1`` would be treated as ``uint16`` data.

        scale_back : bool, optional
            If `True`, when saving changes to a file that contained scaled
            image data, restore the data to the original type and reapply the
            original BSCALE/BZERO values.  This could lead to loss of accuracy
            if scaling back to integer values after performing floating point
            operations on the data.
        """

        # This __init__ currently does nothing differently from the base class,
        # and is only explicitly defined for the docstring.

        super(ImageHDU, self).__init__(
            data=data, header=header, name=name,
            do_not_scale_image_data=do_not_scale_image_data, uint=uint,
            scale_back=scale_back)

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        xtension = card.value
        if isinstance(xtension, string_types):
            xtension = xtension.rstrip()
        return card.keyword == 'XTENSION' and xtension == cls._extension

    def _verify(self, option='warn'):
        """
        ImageHDU verify method.
        """

        errs = super(ImageHDU, self)._verify(option=option)
        naxis = self._header.get('NAXIS', 0)
        # PCOUNT must == 0, GCOUNT must == 1; the former is verified in
        # ExtensionHDU._verify, however ExtensionHDU._verify allows PCOUNT
        # to be >= 0, so we need to check it here
        self.req_cards('PCOUNT', naxis + 3, lambda v: (_is_int(v) and v == 0),
                       0, option, errs)
        return errs


class _IndexInfo(object):
    def __init__(self, indx, naxis):
        if _is_int(indx):
            if 0 <= indx < naxis:
                self.npts = 1
                self.offset = indx
                self.contiguous = True
            else:
                raise IndexError('Index %s out of range.' % indx)
        elif isinstance(indx, slice):
            start, stop, step = indx.indices(naxis)
            self.npts = (stop - start) // step
            self.offset = start
            self.contiguous = step == 1
        elif isiterable(indx):
            self.npts = len(indx)
            self.offset = 0
            self.contiguous = False
        else:
            raise IndexError('Illegal index %s' % indx)
