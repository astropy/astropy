# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division

import collections
import copy
import inspect
import itertools
import os
import re
import sys
import warnings

from .card import Card, CardList, _pad, BLANK_CARD, KEYWORD_LENGTH
from .file import _File, PYTHON_MODES
from .util import (encode_ascii, decode_ascii, fileobj_mode,
                   fileobj_is_binary)

from ...utils import deprecated, isiterable


PY3K = sys.version_info[:2] >= (3, 0)


BLOCK_SIZE = 2880  # the FITS block size

HEADER_END_RE = re.compile(encode_ascii('END {77} *'))


# According to the FITS standard the only characters that may appear in a
# header record are the restricted ASCII chars from 0x20 through 0x7E.
VALID_HEADER_CHARS = set(chr(x) for x in range(0x20, 0x7F))


class Header(object):
    """
    FITS header class.  This class exposes both a dict-like interface and a
    list-like interface to FITS headers.

    The header may be indexed by keyword and, like a dict, the associated value
    will be returned.  When the header contains cards with duplicate keywords,
    only the value of the first card with the given keyword will be returned.
    It is also possible to use a 2-tuple as the index in the form (keyword,
    n)--this returns the n-th value with that keyword, in the case where there
    are duplicate keywords.

    For example:

        >>> header['NAXIS']
        0
        >>> header[('FOO', 1)] # Return the value of the second FOO keyword
        'foo'

    The header may also be indexed by card number:

        >>> header[0] # Return the value of the first card in the header
        'T'

    Commentary keywords such as HISTORY and COMMENT are special cases: When
    indexing the Header object with either 'HISTORY' or 'COMMENT' a list of all
    the HISTORY/COMMENT values is returned:

        >>> header['HISTORY']
        This is the first history entry in this header.
        This is the second history entry in this header.
        ...

    See the Astropy documentation for more details on working with headers.
    """

    def __init__(self, cards=[], txtfile=None):
        """
        Construct a `Header` from an iterable and/or text file.

        Parameters
        ----------
        cards : A list of `Card` objects (optional)
            The cards to initialize the header with.

        txtfile : file path, file object or file-like object (optional)
            Input ASCII header parameters file **(Deprecated)**
            Use the Header.fromfile classmethod instead.
        """

        self.clear()

        if txtfile:
            warnings.warn(
                'The txtfile argument is deprecated.  Use Header.fromfile to '
                'create a new Header object from a text file.',
                DeprecationWarning)
            # get the cards from the input ASCII file
            self.update(self.fromfile(txtfile))
            self._modified = False
            return

        if isinstance(cards, Header):
            cards = cards.cards

        for card in cards:
            self.append(card, end=True)

        self._modified = False

    def __len__(self):
        return len(self._cards)

    def __iter__(self):
        for card in self._cards:
            yield card.keyword

    def __contains__(self, keyword):
        try:
            self._cardindex(keyword)
        except (KeyError, IndexError):
            return False
        return True

    def __getitem__(self, key):
        if isinstance(key, slice):
            return Header([copy.copy(c) for c in self._cards[key]])
        elif self._haswildcard(key):
            return Header([copy.copy(self._cards[idx])
                           for idx in self._wildcardmatch(key)])
        elif (isinstance(key, basestring) and
              key.upper() in Card._commentary_keywords):
            key = key.upper()
            # Special case for commentary cards
            return _HeaderCommentaryCards(self, key)
        if isinstance(key, tuple):
            keyword = key[0]
        else:
            keyword = key
        card = self._cards[self._cardindex(key)]
        if (card.field_specifier is not None and
            keyword == card.keyword.split('.', 1)[0]):
            # This is RVKC; if only the top-level keyword was specified return
            # the raw value, not the parsed out float value
            return card.rawvalue
        return card.value

    def __setitem__(self, key, value):
        if isinstance(key, slice) or self._haswildcard(key):
            if isinstance(key, slice):
                indices = xrange(*key.indices(len(self)))
            else:
                indices = self._wildcardmatch(key)
            if isinstance(value, basestring) or not isiterable(value):
                value = itertools.repeat(value, len(indices))
            for idx, val in itertools.izip(indices, value):
                self[idx] = val
            return

        if isinstance(value, tuple):
            if not (0 < len(value) <= 2):
                raise ValueError(
                    'A Header item may be set with either a scalar value, '
                    'a 1-tuple containing a scalar value, or a 2-tuple '
                    'containing a scalar value and comment string.')
            if len(value) == 1:
                value, comment = value[0], None
                if value is None:
                    value = ''
            elif len(value) == 2:
                value, comment = value
                if value is None:
                    value = ''
                if comment is None:
                    comment = ''
        else:
            comment = None

        card = None
        if isinstance(key, int):
            card = self._cards[key]
        elif isinstance(key, tuple):
            card = self._cards[self._cardindex(key)]
        if card:
            card.value = value
            if comment is not None:
                card.comment = comment
            if card._modified:
                self._modified = True
        else:
            # If we get an IndexError that should be raised; we don't allow
            # assignment to non-existing indices
            self._update((key, value, comment))

    def __delitem__(self, key):
        if isinstance(key, slice) or self._haswildcard(key):
            # This is very inefficient but it's not a commonly used feature.
            # If someone out there complains that they make heavy use of slice
            # deletions and it's too slow, well, we can worry about it then
            # [the solution is not too complicated--it would be wait 'til all
            # the cards are deleted before updating _keyword_indices rather
            # than updating it once for each card that gets deleted]
            if isinstance(key, slice):
                indices = xrange(*key.indices(len(self)))
                # If the slice step is backwards we want to reverse it, because
                # it will be reversed in a few lines...
                if key.step and key.step < 0:
                    indices = reversed(indices)
            else:
                indices = self._wildcardmatch(key)
            for idx in reversed(indices):
                del self[idx]
            return
        elif isinstance(key, basestring):
            # delete ALL cards with the same keyword name
            key = Card.normalize_keyword(key)
            if key not in self._keyword_indices:
                if _is_astropy_internal():
                    # All internal code is designed to assume that this will
                    # raise a KeyError, so go ahead and do so
                    raise KeyError("Keyword '%s' not found." % key)
                # Warn everyone else.
                # TODO: Remove this warning and make KeyError the default after
                # a couple versions (by 3.2 or 3.3, say)
                warnings.warn(
                    'Deletetion of non-existent keyword %r: '
                    'In a future Astropy version Header.__delitem__ may be '
                    'changed so that this raises a KeyError just like a dict '
                    'would. Please update your code so that KeyErrors are '
                    'caught and handled when deleting non-existent keywords.' %
                    key, DeprecationWarning)
                return
            for idx in reversed(self._keyword_indices[key]):
                # Have to copy the indices list since it will be modified below
                del self[idx]
            return

        idx = self._cardindex(key)
        keyword = self._cards[idx].keyword
        del self._cards[idx]
        indices = self._keyword_indices[keyword]
        indices.remove(idx)
        if not indices:
            del self._keyword_indices[keyword]

        # We also need to update all other indices
        self._updateindices(idx, increment=False)
        self._modified = True

    def __repr__(self):
        return self.tostring(sep='\n', endcard=False, padding=False)

    def __str__(self):
        return self.tostring()

    def __eq__(self, other):
        """
        Two Headers are equal only if they have the exact same string
        representation.
        """

        return str(self) == str(other)

    def __add__(self, other):
        temp = self.copy(strip=False)
        temp.extend(other)
        return temp

    def __iadd__(self, other):
        self.extend(other)
        return self

    @property
    def cards(self):
        """
        The underlying physical cards that make up this Header; it can be
        looked at, but it should not be modified directly.
        """

        return _CardAccessor(self)

    @property
    def comments(self):
        """
        View the comments associated with each keyword, if any.

        For example, to see the comment on the NAXIS keyword:

            >>> header.comments['NAXIS']
            number of data axes

        Comments can also be updated through this interface:

            >>> header.comments['NAXIS'] = 'Number of data axes'

        """

        return _HeaderComments(self)

    @property
    def _modified(self):
        """
        Whether or not the header has been modified; this is a property so that
        it can also check each card for modifications--cards may have been
        modified directly without the header containing it otherwise knowing.
        """

        modified_cards = any(c._modified for c in self._cards)
        if modified_cards:
            # If any cards were modified then by definition the header was
            # modified
            self.__dict__['_modified'] = True

        return self.__dict__['_modified']

    @_modified.setter
    def _modified(self, val):
        self.__dict__['_modified'] = val

    @classmethod
    def fromstring(cls, data, sep=''):
        """
        Creates an HDU header from a byte string containing the entire header
        data.

        Parameters
        ----------
        data : str
           String containing the entire header.

        sep : str (optional)
            The string separating cards from each other, such as a newline.  By
            default there is no card separator (as is the case in a raw FITS
            file).

        Returns
        -------
        header
            A new `Header` instance.
        """

        cards = []

        end = 'END' + ' ' * 77

        # If the card separator contains characters that may validly appear in
        # a card, the only way to unambiguously distinguish between cards is to
        # require that they be Card.length long.  However, if the separator
        # contains non-valid characters (namely \n) the cards may be split
        # immediately at the separator
        require_full_cardlength = set(sep).issubset(VALID_HEADER_CHARS)

        # Split the header into individual cards
        idx = 0
        image = []

        while idx < len(data):
            if require_full_cardlength:
                end_idx = idx + Card.length
            else:
                try:
                    end_idx = data.index(sep, idx)
                except ValueError:
                    end_idx = len(data)

            next_image = data[idx:end_idx]
            idx = end_idx + len(sep)

            if image:
                if next_image[:8] == 'CONTINUE':
                    image.append(next_image)
                    continue
                cards.append(Card.fromstring(''.join(image)))

            if require_full_cardlength:
                if next_image == end:
                    image = []
                    break
            else:
                if next_image.split(sep)[0].rstrip() == 'END':
                    image = []
                    break

            image = [next_image]

        # Add the last image that was found before the end, if any
        if image:
            cards.append(Card.fromstring(''.join(image)))

        return cls(cards)

    @classmethod
    def fromfile(cls, fileobj, sep='', endcard=True, padding=True):
        """
        Similar to :meth:`Header.fromstring`, but reads the header string from
        a given file-like object or filename.

        Parameters
        ----------
        fileobj : str, file-like
            A filename or an open file-like object from which a FITS header is
            to be read.  For open file handles the file pointer must be at the
            beginning of the header.

        sep : str (optional)
            The string separating cards from each other, such as a newline.  By
            default there is no card separator (as is the case in a raw FITS
            file).

        endcard : bool (optional)
            If True (the default) the header must end with an END card in order
            to be considered valid.  If an END card is not found an `IOError`
            is raised.

        padding : bool (optional)
            If True (the default) the header will be required to be padded out
            to a multiple of 2880, the FITS header block size.  Otherwise any
            padding, or lack thereof, is ignored.

        Returns
        -------
        header
            A new `Header` instance.
        """

        close_file = False
        if isinstance(fileobj, basestring):
            # Open in text mode by default to support newline handling; if a
            # binary-mode file object is passed in, the user is on their own
            # with respect to newline handling
            fileobj = open(fileobj, 'r')
            close_file = True

        is_binary = fileobj_is_binary(fileobj)
        actual_block_size = _block_size(sep)
        clen = Card.length + len(sep)

        try:
            # Read the first header block.
            block = fileobj.read(actual_block_size)
            if not is_binary:
                block = encode_ascii(block)

            if not block:
                raise EOFError()

            blocks = []
            is_eof = False

            # continue reading header blocks until END card is reached
            while True:
                # find the END card
                is_end = False
                for mo in HEADER_END_RE.finditer(block):
                    # Ensure the END card was found, and it started on the
                    # boundary of a new card (see ticket #142)
                    if mo.start() % clen == 0:
                        # This must be the last header block, otherwise the
                        # file is malformatted
                        is_end = True
                        break

                if not is_end:
                    blocks.append(decode_ascii(block))
                    block = fileobj.read(actual_block_size)
                    if not is_binary:
                        block = encode_ascii(block)
                    if not block:
                        is_eof = True
                        break
                else:
                    break

            last_block = block
            blocks.append(decode_ascii(block))

            blocks = ''.join(blocks)

            # Strip any zero-padding (see ticket #106)
            if blocks and blocks[-1] == '\0':
                if is_eof and blocks.strip('\0') == '':
                    warnings.warn('Unexpected extra padding at the end of the '
                                  'file.  This padding may not be preserved '
                                  'when saving changes.')
                    raise EOFError()
                else:
                    # Replace the illegal null bytes with spaces as required by
                    # the FITS standard, and issue a nasty warning
                    warnings.warn('Header block contains null bytes instead '
                                  'of spaces for padding, and is not FITS-'
                                  'compliant. Nulls may be replaced with '
                                  'spaces upon writing.')
                    blocks.replace('\0', ' ')

            if not HEADER_END_RE.search(last_block) and endcard:
                raise IOError('Header missing END card.')

            if padding and (len(blocks) % actual_block_size) != 0:
                # This error message ignores the length of the separator for
                # now, but maybe it shouldn't?
                actual_len = len(blocks) - actual_block_size + BLOCK_SIZE
                raise ValueError('Header size is not multiple of %d: %d'
                                 % (BLOCK_SIZE, actual_len))

            return cls.fromstring(blocks, sep=sep)
        finally:
            if close_file:
                fileobj.close()

    def tostring(self, sep='', endcard=True, padding=True):
        r"""
        Returns a string representation of the header.

        By default this uses no separator between cards, adds the END card, and
        pads the string with spaces to the next multiple of 2880 bytes.  That
        is, it returns the header exactly as it would appear in a FITS file.

        Parameters
        ----------
        sep : str (optional)
            The character or string with which to separate cards.  By default
            there is no separator, but one could use `'\\n'`, for example, to
            separate each card with a new line

        endcard : bool (optional)
            If True (default) adds the END card to the end of the header
            string

        padding : bool (optional)
            If True (default) pads the string with spaces out to the next
            multiple of 2880 characters

        Returns
        -------
        s : string
            A string representing a FITS header.
        """

        lines = []
        for card in self._cards:
            s = str(card)
            # Cards with CONTINUE cards may be longer than 80 chars; so break
            # them into multiple lines
            while s:
                lines.append(s[:Card.length])
                s = s[Card.length:]

        s = sep.join(lines)
        if endcard:
            s += sep + _pad('END')
        if padding:
            s += ' ' * _pad_length(len(s))
        return s

    def tofile(self, fileobj, sep='', endcard=True, padding=True,
               clobber=False):
        r"""
        Writes the header to file or file-like object.

        By default this writes the header exactly as it would be written to a
        FITS file, with the END card included and padding to the next multiple
        of 2880 bytes.  However, aspects of this may be controlled.

        Parameters
        ----------
        fileobj : str, file (optional)
            Either the pathname of a file, or an open file handle or file-like
            object

        sep : str (optional)
            The character or string with which to separate cards.  By default
            there is no separator, but one could use `'\\n'`, for example, to
            separate each card with a new line

        endcard : bool (optional)
            If `True` (default) adds the END card to the end of the header
            string

        padding : bool (optional)
            If `True` (default) pads the string with spaces out to the next
            multiple of 2880 characters

        clobber : bool (optional)
            If `True`, overwrites the output file if it already exists
        """

        close_file = False

        # check if the output file already exists
        # TODO: Perhaps this sort of thing could be handled by the _File
        # initializer...
        if isinstance(fileobj, basestring):
            if os.path.exists(fileobj) and os.path.getsize(fileobj) != 0:
                if clobber:
                    warnings.warn("Overwriting existing file '%s'." % fileobj)
                    os.remove(fileobj)
                else:
                    raise IOError("File '%s' already exists." % fileobj)

            fileobj = open(fileobj, 'wb')
            close_file = True

        if not isinstance(fileobj, _File):
            # TODO: There needs to be a way of handling this built into the
            # _File class.  I think maybe there used to be, but I took it out;
            # now the design is such that it would be better for it to go back
            # in
            mode = 'append'
            fmode = fileobj_mode(fileobj) or 'ab+'
            for key, val in PYTHON_MODES.iteritems():
                if val == fmode:
                    mode = key
                    break
            fileobj = _File(fileobj, mode=mode)

        try:
            blocks = self.tostring(sep=sep, endcard=endcard, padding=padding)
            actual_block_size = _block_size(sep)
            if padding and len(blocks) % actual_block_size != 0:
                raise IOError('Header size (%d) is not a multiple of block '
                              'size (%d).' %
                              (len(blocks) - actual_block_size + BLOCK_SIZE,
                               BLOCK_SIZE))

            if not fileobj.simulateonly:
                fileobj.flush()
                try:
                    offset = fileobj.tell()
                except (AttributeError, IOError):
                    offset = 0
                fileobj.write(blocks.encode('ascii'))
                fileobj.flush()
        finally:
            if close_file:
                fileobj.close()

    @classmethod
    def fromtextfile(cls, fileobj, endcard=False):
        """
        Equivalent to ``Header.fromfile(fileobj, sep='\\n', endcard=False,
        padding=False)``.
        """

        return cls.fromfile(fileobj, sep='\n', endcard=endcard, padding=False)

    def totextfile(self, fileobj, endcard=False, clobber=False):
        """
        Equivalent to ``Header.tofile(fileobj, sep='\\n', endcard=False,
        padding=False, clobber=clobber)``.
        """

        self.tofile(fileobj, sep='\n', endcard=endcard, padding=False,
                    clobber=clobber)

    def clear(self):
        """
        Remove all cards from the header.
        """

        self._cards = []
        self._keyword_indices = collections.defaultdict(list)

    def copy(self, strip=False):
        """
        Make a copy of the :class:`Header`.

        Parameters
        ----------
        strip : bool (optional)
           If True, strip any headers that are specific to one of the standard
           HDU types, so that this header can be used in a different HDU.

        Returns
        -------
        header
            A new :class:`Header` instance.
        """

        tmp = Header([copy.copy(card) for card in self._cards])
        if strip:
            tmp._strip()
        return tmp

    @classmethod
    def fromkeys(cls, iterable, value=None):
        """
        Similar to :meth:`dict.fromkeys`--creates a new `Header` from an
        iterable of keywords and an optional default value.

        This method is not likely to be particularly useful for creating real
        world FITS headers, but it is useful for testing.

        Parameters
        ----------
        iterable
            Any iterable that returns strings representing FITS keywords.

        value : (optional)
            A default value to assign to each keyword; must be a valid type for
            FITS keywords.

        Returns
        -------
        header
            A new `Header` instance.
        """

        d = cls()
        if not isinstance(value, tuple):
            value = (value,)
        for key in iterable:
            d.append((key,) + value)
        return d

    def get(self, key, default=None):
        """
        Similar to :meth:`dict.get`--returns the value associated with keyword
        in the header, or a default value if the keyword is not found.

        Parameters
        ----------
        key : str
            A keyword that may or may not be in the header.

        default : (optional)
            A default value to return if the keyword is not found in the
            header.

        Returns
        -------
        value
            The value associated with the given keyword, or the default value
            if the keyword is not in the header.
        """

        try:
            return self[key]
        except (KeyError, IndexError):
            return default

    def set(self, keyword, value=None, comment=None, before=None, after=None):
        """
        Set the value and/or comment and/or position of a specified keyword.

        If the keyword does not already exist in the header, a new keyword is
        created in the specified position, or appended to the end of the header
        if no position is specified.

        This method is similar to :meth:`Header.update` prior to PyFITS 3.1.

        .. note::
            It should be noted that ``header.set(keyword, value)`` and
            ``header.set(keyword, value, comment)`` are equivalent to
            ``header[keyword] = value`` and
            ``header[keyword] = (value, comment)`` respectfully.

            New keywords can also be inserted relative to existing keywords
            using, for example
            ``header.insert('NAXIS1', ('NAXIS', 2, 'Number of axes'))`` to
            insert before an existing keyword, or
            ``header.insert('NAXIS', ('NAXIS1', 4096), after=True)`` to insert
            after an existing keyword.

            The the only advantage of using :meth:`Header.set` is that it
            easily replaces the old usage of :meth:`Header.update` both
            conceptually and in terms of function signature.

        Parameters
        ----------
        keyword : str
            A header keyword

        value : str (optional)
            The value to set for the given keyword; if None the existing value
            is kept, but '' may be used to set a blank value

        comment : str (optional)
            The comment to set for the given keyword; if None the existing
            comment is kept, but '' may be used to set a blank comment

        before : str, int (optional)
            Name of the keyword, or index of the `Card` before which
            this card should be located in the header.  The argument `before`
            takes precedence over `after` if both specified.

        after : str, int (optional)
            Name of the keyword, or index of the `Card` after which this card
            should be located in the header.

        """

        # Create a temporary card that looks like the one being set; if the
        # temporary card turns out to be a RVKC this will make it easier to
        # deal with the idiosyncrasies thereof
        # Don't try to make a temporary card though if they keyword looks like
        # it might be a HIERARCH card or is otherwise invalid--this step is
        # only for validating RVKCs.
        if (len(keyword) <= KEYWORD_LENGTH and
            Card._keywd_FSC_RE.match(keyword) and
            keyword not in self._keyword_indices):
            new_card = Card(keyword, value, comment)
            new_keyword = new_card.keyword
        else:
            new_keyword = keyword

        if (new_keyword not in Card._commentary_keywords and
            new_keyword in self):
            if comment is None:
                comment = self.comments[keyword]
            if value is None:
                value = self[keyword]

            self[keyword] = (value, comment)

            if before is not None or after is not None:
                card = self._cards[self._cardindex(keyword)]
                self._relativeinsert(card, before=before, after=after,
                                     replace=True)
        elif before is not None or after is not None:
            self._relativeinsert((keyword, value, comment), before=before,
                                 after=after)
        else:
            self[keyword] = (value, comment)

    @deprecated('3.0', alternative='``key in header`` syntax')
    def has_key(self, key):
        """Like :meth:`dict.has_key`."""

        return key in self

    def items(self):
        """Like :meth:`dict.items`."""

        return list(self.iteritems())

    def iteritems(self):
        """Like :meth:`dict.iteritems`."""

        for card in self._cards:
            yield (card.keyword, card.value)

    def iterkeys(self):
        """
        Like :meth:`dict.iterkeys`--iterating directly over the `Header`
        instance has the same behavior.
        """

        return self.__iter__()

    def itervalues(self):
        """Like :meth:`dict.itervalues`."""

        for _, v in self.iteritems():
            yield v

    def keys(self):
        """
        Return a list of keywords in the header in the order they
        appear--like:meth:`dict.keys` but ordered.
        """

        return [keyword for keyword in self]

    def pop(self, *args):
        """
        Works like :meth:`list.pop` if no arguments or an index argument are
        supplied; otherwise works like :meth:`dict.pop`.
        """

        if len(args) > 2:
            raise TypeError('Header.pop expected at most 2 arguments, got '
                            '%d' % len(args))

        if len(args) == 0:
            key = -1
        else:
            key = args[0]

        try:
            value = self[key]
        except (KeyError, IndexError):
            if len(args) == 2:
                return args[1]
            raise

        del self[key]
        return value

    def popitem(self):
        try:
            k, v = self.iteritems().next()
        except StopIteration:
            raise KeyError('Header is empty')
        del self[k]
        return k, v

    def setdefault(self, key, default=None):
        try:
            return self[key]
        except (KeyError, IndexError):
            self[key] = default
        return default

    def update(self, *args, **kwargs):
        """
        Update the Header with new keyword values, updating the values of
        existing keywords and appending new keywords otherwise; similar to
        dict.update().

        update() accepts either a dict-like object or an iterable.  In the
        former case the keys must be header keywords and the values may be
        either scalar values or (value, comment) tuples.  In the case of an
        iterable the items must be (keyword, value) tuples or
        (keyword, value, comment) tuples.

        Arbitrary arguments are also accepted, in which case the update() is
        called again with the kwargs dict as its only argument.  That is,

            >>> header.update(NAXIS1=100, NAXIS2=100)

        is equivalent to

            >>> header.update({'NAXIS1': 100, 'NAXIS2': 100})

        .. warning::
            As this method works similarly to dict.update() it is very
            different from the Header.update() method in PyFITS versions prior
            to 3.1.0.  However, support for the old API is also maintained for
            backwards compatibility.  If update() is called with at least two
            positional arguments then it can be assumed that the old API is
            being used.  Use of the old API should be considered
            **deprecated**.  Most uses of the old API can be replaced as
            follows:

            * Replace

                  >>> header.update(keyword, value)

              with

                  >>> header[keyword] = value

            * Replace

                  >>> header.update(keyword, value, comment=comment)

              with

                  >>> header[keyword] = (value, comment)

            * Replace

                  >>> header.update(keyword, value, before=before_keyword)

              with

                  >>> header.insert(before_keyword, (keyword, value))

            * Replace

                  >>> header.update(keyword, value, after=after_keyword)

              with

                  >>> header.insert(after_keyword, (keyword, value),
                  ...               after=True)

            See also :meth:`Header.set` which is a new method that provides an
            interface similar to the old Header.update() and may help make
            transition a little easier.

            For reference, the old documentation for the old Header.update()
            is provided below:

        Update one header card.

        If the keyword already exists, it's value and/or comment will
        be updated.  If it does not exist, a new card will be created
        and it will be placed before or after the specified location.
        If no `before` or `after` is specified, it will be appended at
        the end.

        Parameters
        ----------
        key : str
            keyword

        value : str
            value to be used for updating

        comment : str (optional)
            to be used for updating, default=None.

        before : str, int (optional)
            name of the keyword, or index of the `Card` before which
            the new card will be placed.  The argument `before` takes
            precedence over `after` if both specified.

        after : str, int (optional)
            name of the keyword, or index of the `Card` after which
            the new card will be placed.

        savecomment : bool (optional)
            When `True`, preserve the current comment for an existing
            keyword.  The argument `savecomment` takes precedence over
            `comment` if both specified.  If `comment` is not
            specified then the current comment will automatically be
            preserved.

        """

        legacy_args = ['key', 'value', 'comment', 'before', 'after',
                       'savecomment']

        # This if statement covers all the cases in which this could be a
        # legacy update(); note that it means it's impossible to do a
        # dict-style update where *all* the keywords happen to legacy
        # arguments, but realistically speaking that use case will not come up

        # The fact that Python is "flexible" in allowing positional args to be
        # passed in as keyword args makes this a little more complicated than
        # it otherwise would be :/
        issubset = set(kwargs).issubset(set(legacy_args))
        if (len(args) >= 2 or
            (len(args) == 1 and 'value' in kwargs and issubset) or
            (len(args) == 0 and 'key' in kwargs and 'value' in kwargs and
             issubset)):
            # This must be a legacy update()
            warnings.warn(
                "The use of header.update() to add new keywords to a header "
                "deprecated.  Instead, use either header.set() or simply "
                "`header[keyword] = value` or "
                "`header[keyword] = (value, comment)`.  header.set() is only "
                "necessary to use if you also want to use the before/after "
                "keyword arguments.", DeprecationWarning)

            for k, v in zip(legacy_args, args):
                if k in kwargs:
                    raise TypeError(
                        '%s.update() got multiple values for keyword '
                        'argument %r' % (self.__class__.__name__, k))
                kwargs[k] = v

            keyword = kwargs.get('key')
            value = kwargs.get('value')
            comment = kwargs.get('comment')
            before = kwargs.get('before')
            after = kwargs.get('after')
            savecomment = kwargs.get('savecomment')

            # Handle the savecomment argument which is not currently used by
            # Header.set()
            if keyword in self and savecomment:
                comment = None

            self.set(keyword, value, comment, before, after)
        else:
            # The rest of this should work similarly to dict.update()
            if args:
                other = args[0]
            else:
                other = None

            def update_from_dict(k, v):
                if not isinstance(v, tuple):
                    card = Card(k, v)
                elif 0 < len(v) <= 2:
                    card = Card(*((k,) + v))
                else:
                    raise ValueError(
                            'Header update value for key %r is invalid; the '
                            'value must be either a scalar, a 1-tuple '
                            'containing the scalar value, or a 2-tuple '
                            'containing the value and a comment string.' % k)
                self._update(card)

            if other is None:
                pass
            elif hasattr(other, 'iteritems'):
                for k, v in other.iteritems():
                    update_from_dict(k, v)
            elif hasattr(other, 'keys'):
                for k in other.keys():
                    update_from_dict(k, other[k])
            else:
                for idx, card in enumerate(other):
                    if isinstance(card, Card):
                        self._update(card)
                    elif isinstance(card, tuple) and (1 < len(card) <= 3):
                        self._update(Card(*card))
                    else:
                        raise ValueError(
                                'Header update sequence item #%d is invalid; '
                                'the item must either be a 2-tuple containing '
                                'a keyword and value, or a 3-tuple containing '
                                'a keyword, value, and comment string.' % idx)
            if kwargs:
                self.update(kwargs)

    def values(self):
        """Returns a list of the values of all cards in the header."""

        return [v for _, v in self.iteritems()]

    def append(self, card=None, useblanks=True, bottom=False, end=False):
        """
        Appends a new keyword+value card to the end of the Header, similar
        to list.append().

        By default if the last cards in the Header have commentary keywords,
        this will append the new keyword before the commentary (unless the new
        keyword is also commentary).

        Also differs from list.append() in that it can be called with no
        arguments: In this case a blank card is appended to the end of the
        Header.  In the case all the keyword arguments are ignored.

        Parameters
        ----------
        card : str, tuple
            A keyword or a (keyword, value, [comment]) tuple representing a
            single header card; the comment is optional in which case a
            2-tuple may be used

        useblanks : bool (optional)
            If there are blank cards at the end of the Header, replace the
            first blank card so that the total number of cards in the Header
            does not increase.  Otherwise preserve the number of blank cards.

        bottom : bool (optional)
            If True, instead of appending after the last non-commentary card,
            append after the last non-blank card.

        end : bool (optional):
            If True, ignore the useblanks and bottom options, and append at the
            very end of the Header.

        """

        if isinstance(card, basestring):
            card = Card(card)
        elif isinstance(card, tuple):
            card = Card(*card)
        elif card is None:
            card = Card()
        elif not isinstance(card, Card):
            raise ValueError(
                'The value appended to a Header must be either a keyword or '
                '(keyword, value, [comment]) tuple; got: %r' % card)

        if not end and str(card) == BLANK_CARD:
            # Blank cards should always just be appended to the end
            end = True

        if end:
            self._cards.append(card)
            idx = len(self._cards) - 1
        else:
            idx = len(self._cards) - 1
            while idx >= 0 and str(self._cards[idx]) == BLANK_CARD:
                idx -= 1

            if not bottom and card.keyword not in Card._commentary_keywords:
                while (idx >= 0 and
                       self._cards[idx].keyword in Card._commentary_keywords):
                    idx -= 1

            idx += 1
            self._cards.insert(idx, card)
            self._updateindices(idx)

        keyword = Card.normalize_keyword(card.keyword)
        self._keyword_indices[keyword].append(idx)

        if not end:
            # If the appended card was a commentary card, and it was appended
            # before existing cards with the same keyword, the indices for
            # cards with that keyword may have changed
            if not bottom and card.keyword in Card._commentary_keywords:
                self._keyword_indices[keyword].sort()

            # Finally, if useblanks, delete a blank cards from the end
            if useblanks:
                self._useblanks(len(str(card)) // Card.length)

        self._modified = True

    def extend(self, cards, strip=True, unique=False, update=False,
               update_first=False, useblanks=True, bottom=False, end=False):
        """
        Appends multiple keyword+value cards to the end of the header, similar
        to list.extend().

        Parameters
        ----------
        cards : iterable
            An iterable of (keyword, value, [comment]) tuples; see
            Header.append()

        strip : bool (optional)
            Remove any keywords that have meaning only to specific types of
            HDUs, so that only more general keywords are added from extension
            Header or Card list (default: True).

        unique : bool (optional)
            If `True`, ensures that no duplicate keywords are appended;
            keywords already in this header are simply discarded.  The
            exception is commentary keywords (COMMENT, HISTORY, etc.): they are
            only treated as duplicates if their values match.

        update : bool (optional)
            If `True`, update the current header with the values and comments
            from duplicate keywords in the input header.  This supercedes the
            `unique` argument.  Commentary keywords are treated the same as if
            `unique=True`.

        update_first : bool (optional)
            If the first keyword in the header is 'SIMPLE', and the first
            keyword in the input header is 'XTENSION', the 'SIMPLE' keyword is
            replaced by the 'XTENSION' keyword.  Likewise if the first keyword
            in the header is 'XTENSION' and the first keyword in the input
            header is 'SIMPLE', the 'XTENSION' keyword is replaced by the
            'SIMPLE' keyword.  This behavior is otherwise dumb as to whether or
            not the resulting header is a valid primary or extension header.
            This is mostly provided to support backwards compatibility with the
            old :meth:`Header.fromTxtFile` method, and only applies if
            `update=True`.

        useblanks, bottom, end : bool (optional)
            These arguments are passed to :meth:`Header.append` while appending
            new cards to the header.
        """

        temp = Header(cards)
        if strip:
            temp._strip()

        if len(self):
            first = self.cards[0].keyword
        else:
            first = None

        # This copy is used to check for duplicates in this header prior to the
        # extend, while not counting duplicates in the header being extended
        # from (see ticket #156)
        orig = self[:]

        for idx, card in enumerate(temp.cards):
            keyword = card.keyword
            if keyword not in Card._commentary_keywords:
                if unique and not update and keyword in orig:
                    continue
                elif update:
                    if idx == 0 and update_first:
                        # Dumbly update the first keyword to either SIMPLE or
                        # XTENSION as the case may be, as was in the case in
                        # Header.fromTxtFile
                        if ((keyword == 'SIMPLE' and first == 'XTENSION') or
                            (keyword == 'XTENSION' and first == 'SIMPLE')):
                            del self[0]
                            self.insert(0, card)
                        else:
                            self[keyword] = (card.value, card.comment)
                    elif keyword in orig:
                        self[keyword] = (card.value, card.comment)
                    else:
                        self.append(card, useblanks=useblanks, bottom=bottom,
                                    end=end)
                else:
                    self.append(card, useblanks=useblanks, bottom=bottom,
                                end=end)
            else:
                if unique or update and keyword in orig:
                    if str(card) == BLANK_CARD:
                        self.append(card, useblanks=useblanks, bottom=bottom,
                                    end=end)
                        continue

                    for value in orig[keyword]:
                        if value == card.value:
                            break
                    else:
                        self.append(card, useblanks=useblanks, bottom=bottom,
                                    end=end)
                else:
                    self.append(card, useblanks=useblanks, bottom=bottom,
                                end=end)

    def count(self, keyword):
        """
        Returns the count of the given keyword in the header, similar to
        list.count() if the Header object is treated as a list of keywords.

        Parameters
        ----------
        keyword : str
            The keyword to count instances of in the header

        """

        keyword = Card.normalize_keyword(keyword)

        # We have to look before we leap, since otherwise _keyword_indices,
        # being a defaultdict, will create an entry for the nonexistent keyword
        if keyword not in self._keyword_indices:
            raise KeyError("Keyword %r not found." % keyword)

        return len(self._keyword_indices[keyword])

    def index(self, keyword, start=None, stop=None):
        """
        Returns the index if the first instance of the given keyword in the
        header, similar to list.index() if the Header object is treated as a
        list of keywords.

        Parameters
        ----------
        keyword : str
            The keyword to look up in the list of all keywords in the header

        start : int (optional)
            The lower bound for the index

        stop : int (optional)
            The upper bound for the index

        """

        if start is None:
            start = 0

        if stop is None:
            stop = len(self._cards)

        if stop < start:
            step = -1
        else:
            step = 1

        keyword = Card.normalize_keyword(keyword)

        for idx in xrange(start, stop, step):
            if self._cards[idx].keyword == keyword:
                return idx
        else:
            raise ValueError('The keyword %r is not in the header.' % keyword)

    def insert(self, idx, card, useblanks=True):
        """
        Inserts a new keyword+value card into the Header at a given location,
        similar to list.insert().

        Parameters
        ----------
        idx : int
            The index into the the list of header keywords before which the
            new keyword should be inserted

        card : str, tuple
            A keyword or a (keyword, value, [comment]) tuple; see
            Header.append()

        useblanks : bool (optional)
            If there are blank cards at the end of the Header, replace the
            first blank card so that the total number of cards in the Header
            does not increase.  Otherwise preserve the number of blank cards.

        """

        if idx >= len(self._cards):
            # This is just an append (Though it must be an append absolutely to
            # the bottom, ignoring blanks, etc.--the point of the insert method
            # is that you get exactly what you asked for with no surprises)
            self.append(card, end=True)
            return

        if isinstance(card, basestring):
            card = Card(card)
        elif isinstance(card, tuple):
            card = Card(*card)
        elif not isinstance(card, Card):
            raise ValueError(
                'The value inserted into a Header must be either a keyword or '
                '(keyword, value, [comment]) tuple; got: %r' % card)

        self._cards.insert(idx, card)

        keyword = card.keyword

        # If idx was < 0, determine the actual index according to the rules
        # used by list.insert()
        if idx < 0:
            idx += len(self._cards) - 1
            if idx < 0:
                idx = 0

        # All the keyword indices above the insertion point must be updated
        self._updateindices(idx)

        self._keyword_indices[keyword].append(idx)
        count = len(self._keyword_indices[keyword])
        if count > 1:
            # There were already keywords with this same name
            if keyword not in Card._commentary_keywords:
                warnings.warn(
                    'A %r keyword already exists in this header.  Inserting '
                    'duplicate keyword.' % keyword)
            self._keyword_indices[keyword].sort()

        if useblanks:
            self._useblanks(len(str(card)) // Card.length)

        self._modified = True

    def remove(self, keyword):
        """
        Removes the first instance of the given keyword from the header
        similar to list.remove() if the Header object is treated as a list of
        keywords.

        Parameters
        ----------
        value : str
            The keyword of which to remove the first instance in the header

        """

        del self[self.index(keyword)]

    def rename_keyword(self, oldkeyword, newkeyword, force=False):
        """
        Rename a card's keyword in the header.

        Parameters
        ----------
        oldkeyword : str or int
            Old keyword or card index

        newkeyword : str
            New keyword

        force : bool (optional)
            When `True`, if the new keyword already exists in the header, force
            the creation of a duplicate keyword.  Otherwise a `ValueError` is
            raised.
        """

        oldkeyword = Card.normalize_keyword(oldkeyword)
        newkeyword = Card.normalize_keyword(newkeyword)

        if newkeyword == 'CONTINUE':
            raise ValueError('Can not rename to CONTINUE')

        if (newkeyword in Card._commentary_keywords or
            oldkeyword in Card._commentary_keywords):
            if not (newkeyword in Card._commentary_keywords and
                    oldkeyword in Card._commentary_keywords):
                raise ValueError('Regular and commentary keys can not be '
                                 'renamed to each other.')
        elif not force and newkeyword in self:
            raise ValueError('Intended keyword %s already exists in header.'
                             % newkeyword)

        idx = self.index(oldkeyword)
        card = self.cards[idx]
        del self[idx]
        self.insert(idx, (newkeyword, card.value, card.comment))

    def add_history(self, value, before=None, after=None):
        """
        Add a ``HISTORY`` card.

        Parameters
        ----------
        value : str
            history text to be added.

        before : str or int, optional
            same as in `Header.update`

        after : str or int, optional
            same as in `Header.update`
        """

        self._add_commentary('HISTORY', value, before=before, after=after)

    def add_comment(self, value, before=None, after=None):
        """
        Add a ``COMMENT`` card.

        Parameters
        ----------
        value : str
            text to be added.

        before : str or int, optional
            same as in `Header.update`

        after : str or int, optional
            same as in `Header.update`
        """

        self._add_commentary('COMMENT', value, before=before, after=after)

    def add_blank(self, value='', before=None, after=None):
        """
        Add a blank card.

        Parameters
        ----------
        value : str, optional
            text to be added.

        before : str or int, optional
            same as in `Header.update`

        after : str or int, optional
            same as in `Header.update`
        """

        self._add_commentary('', value, before=before, after=after)

    def _update(self, card):
        """
        The real update code.  If keyword already exists, its value and/or
        comment will be updated.  Otherwise a new card will be appended.

        This will not create a duplicate keyword except in the case of
        commentary cards.  The only other way to force creation of a duplicate
        is to use the insert(), append(), or extend() methods.
        """

        keyword, value, comment = card

        # Lookups for existing/known keywords are case-insensitive
        keyword = keyword.upper()
        if keyword.startswith('HIERARCH '):
            keyword = keyword[9:]

        if (keyword not in Card._commentary_keywords and
            keyword in self._keyword_indices):
            # Easy; just update the value/comment
            idx = self._keyword_indices[keyword][0]
            existing_card = self._cards[idx]
            existing_card.value = value
            if comment is not None:
                # '' should be used to explictly blank a comment
                existing_card.comment = comment
            if existing_card._modified:
                self._modified = True
        elif keyword in Card._commentary_keywords:
            cards = self._splitcommentary(keyword, value)
            if keyword in self._keyword_indices:
                # Append after the last keyword of the same type
                idx = self.index(keyword, start=len(self) - 1, stop=-1)
                isblank = not (keyword or value or comment)
                for c in reversed(cards):
                    self.insert(idx + 1, c, useblanks=(not isblank))
            else:
                for c in cards:
                    self.append(c, bottom=True)
        else:
            # A new keyword! self.append() will handle updating _modified
            self.append(card)

    def _cardindex(self, key):
        """Returns an index into the ._cards list given a valid lookup key."""

        if isinstance(key, slice):
            return key
        elif isinstance(key, int):
            # If < 0, determine the actual index
            if key < 0:
                key += len(self._cards)
            if key < 0 or key >= len(self._cards):
                raise IndexError('Header index out of range.')
            return key

        if isinstance(key, basestring):
            key = (key, 0)

        if isinstance(key, tuple):
            if (len(key) != 2 or not isinstance(key[0], basestring) or
                    not isinstance(key[1], int)):
                raise ValueError(
                        'Tuple indices must be 2-tuples consisting of a '
                        'keyword string and an integer index.')
            keyword, n = key
            keyword = Card.normalize_keyword(keyword)
            # Returns the index into _cards for the n-th card with the given
            # keyword (where n is 0-based)
            if keyword and keyword not in self._keyword_indices:
                if len(keyword) > KEYWORD_LENGTH or '.' in keyword:
                    raise KeyError("Keyword %r not found." % keyword)
                # Great--now we have to check if there's a RVKC that starts
                # with the given keyword, making failed lookups fairly
                # expensive
                # TODO: Find a way to make this more efficient; perhaps a set
                # of RVKCs in the header or somesuch.
                keyword = keyword + '.'
                found = 0
                for idx, card in enumerate(self._cards):
                    if (card.field_specifier and
                        card.keyword.startswith(keyword)):
                        if found == n:
                            return idx
                        found += 1
                else:
                    raise KeyError("Keyword %r not found." % keyword[:-1])
            try:
                return self._keyword_indices[keyword][n]
            except IndexError:
                raise IndexError('There are only %d %r cards in the header.' %
                                 (len(self._keyword_indices[keyword]),
                                  keyword))
        else:
            raise ValueError(
                    'Header indices must be either a string, a 2-tuple, or '
                    'an integer.')

    def _relativeinsert(self, card, before=None, after=None, replace=False):
        """
        Inserts a new card before or after an existing card; used to
        implement support for the legacy before/after keyword arguments to
        Header.update().

        If replace=True, move an existing card with the same keyword.
        """

        if before is None:
            insertionkey = after
        else:
            insertionkey = before

        def get_insertion_idx():
            if not (isinstance(insertionkey, int) and
                    insertionkey >= len(self._cards)):
                idx = self._cardindex(insertionkey)
            else:
                idx = insertionkey

            if before is None:
                idx += 1

            return idx

        if replace:
            # The card presumably already exists somewhere in the header.
            # Check whether or not we actually have to move it; if it does need
            # to be moved we just delete it and then it will be reinstered
            # below
            old_idx = self._cardindex(card.keyword)
            insertion_idx = get_insertion_idx()

            if (insertion_idx >= len(self._cards) and
                old_idx == len(self._cards) - 1):
                # The card would be appended to the end, but it's already at
                # the end
                return

            if before is not None:
                if old_idx == insertion_idx - 1:
                    return
            elif after is not None and old_idx == insertion_idx:
                return

            del self[old_idx]


        # Even if replace=True, the insertion idx may have changed since the
        # old card was deleted
        idx = get_insertion_idx()

        if card[0] in Card._commentary_keywords:
            cards = reversed(self._splitcommentary(card[0], card[1]))
        else:
            cards = [card]
        for c in cards:
            self.insert(idx, c)

    def _updateindices(self, idx, increment=True):
        """
        For all cards with index above idx, increment or decrement its index
        value in the keyword_indices dict.
        """

        if idx > len(self._cards):
            # Save us some effort
            return

        increment = 1 if increment else -1

        for indices in self._keyword_indices.itervalues():
            for jdx, keyword_index in enumerate(indices):
                if keyword_index >= idx:
                    indices[jdx] += increment

    def _countblanks(self):
        """Returns the number of blank cards at the end of the Header."""

        for idx in xrange(1, len(self._cards)):
            if str(self._cards[-idx]) != BLANK_CARD:
                return idx - 1
        return 0

    def _useblanks(self, count):
        for _ in range(count):
            if str(self._cards[-1]) == BLANK_CARD:
                del self[-1]
            else:
                break

    def _haswildcard(self, keyword):
        """Return `True` if the input keyword contains a wildcard pattern."""

        return (isinstance(keyword, basestring) and
                (keyword.endswith('...') or '*' in keyword or '?' in keyword))

    def _wildcardmatch(self, pattern):
        """
        Returns a list of indices of the cards matching the given wildcard
        pattern.

         * '*' matches 0 or more alphanumeric characters or _
         * '?' matches a single alphanumeric character or _
         * '...' matches 0 or more of any non-whitespace character
        """

        pattern = pattern.replace('*', r'\w*').replace('?', r'\w')
        pattern = pattern.replace('...', r'\S*') + '$'
        pattern_re = re.compile(pattern, re.I)

        return [idx for idx, card in enumerate(self._cards)
                if pattern_re.match(card.keyword)]

    def _splitcommentary(self, keyword, value):
        """
        Given a commentary keyword and value, returns a list of the one or more
        cards needed to represent the full value.  This is primarily used to
        create the multiple commentary cards needed to represent a long value
        that won't fit into a single commentary card.
        """

        # The maximum value in each card can be the maximum card length minus
        # the maximum key length (which can include spaces if they key length
        # less than 8
        maxlen = Card.length - KEYWORD_LENGTH
        valuestr = str(value)

        if len(valuestr) <= maxlen:
            # The value can fit in a single card
            cards = [Card(keyword, value)]
        else:
            # The value must be split across multiple consecutive commentary
            # cards
            idx = 0
            cards = []
            while idx < len(valuestr):
                cards.append(Card(keyword, valuestr[idx:idx + maxlen]))
                idx += maxlen
        return cards

    def _strip(self):
        """
        Strip cards specific to a certain kind of header.

        Strip cards like ``SIMPLE``, ``BITPIX``, etc. so the rest of
        the header can be used to reconstruct another kind of header.
        """

        # TODO: Previously this only deleted some cards specific to an HDU if
        # _hdutype matched that type.  But it seemed simple enough to just
        # delete all desired cards anyways, and just ignore the KeyErrors if
        # they don't exist.
        # However, it might be desirable to make this extendable somehow--have
        # a way for HDU classes to specify some headers that are specific only
        # to that type, and should be removed otherwise.

        if 'NAXIS' in self:
            naxis = self['NAXIS']
        else:
            naxis = 0

        if 'TFIELDS' in self:
            tfields = self['TFIELDS']
        else:
            tfields = 0

        for idx in range(naxis):
            try:
                del self['NAXIS' + str(idx + 1)]
            except KeyError:
                pass

        for name in ('TFORM', 'TSCAL', 'TZERO', 'TNULL', 'TTYPE',
                     'TUNIT', 'TDISP', 'TDIM', 'THEAP', 'TBCOL'):
            for idx in range(tfields):
                try:
                    del self[name + str(idx + 1)]
                except KeyError:
                    pass

        for name in ('SIMPLE', 'XTENSION', 'BITPIX', 'NAXIS', 'EXTEND',
                     'PCOUNT', 'GCOUNT', 'GROUPS', 'BSCALE', 'BZERO',
                     'TFIELDS'):
            try:
                del self[name]
            except KeyError:
                pass

    def _add_commentary(self, key, value, before=None, after=None):
        """
        Add a commentary card.

        If `before` and `after` are `None`, add to the last occurrence
        of cards of the same name (except blank card).  If there is no
        card (or blank card), append at the end.
        """

        if before is not None or after is not None:
            self._relativeinsert((key, value), before=before,
                                 after=after)
        else:
            self[key] = value

    # Some fixes for compatibility with the Python 3 dict interface, where
    # iteritems -> items, etc.
    if PY3K:  # pragma: py3
        keys = iterkeys
        values = itervalues
        items = iteritems
        del iterkeys
        del itervalues
        del iteritems

    # The following properties/methods are for legacy API backwards
    # compatibility

    @property
    @deprecated('3.1', alternative='the `.cards` attribute')
    def ascard(self):
        """
        Returns a CardList object wrapping this Header; provided for
        backwards compatibility for the old API (where Headers had an
        underlying CardList).
        """

        return CardList(self)

    @deprecated('3.0', alternative='the `.ascard` attribute')
    def ascardlist(self):
        """
        Returns a `CardList` object.
        """

        return self.ascard

    @deprecated('3.1', alternative=':meth:`Header.rename_keyword`')
    def rename_key(self, oldkey, newkey, force=False):
        self.rename_keyword(oldkey, newkey, force)

    @deprecated('3.1', alternative="``header['HISTORY']``", pending=True)
    def get_history(self):
        """
        Get all history cards as a list of string texts.
        """

        if 'HISTORY' in self:
            return self['HISTORY']
        else:
            return []

    @deprecated('3.1', alternative="``header['COMMENT']``", pending=True)
    def get_comment(self):
        """
        Get all comment cards as a list of string texts.
        """

        if 'COMMENT' in self:
            return self['COMMENT']
        else:
            return []

    @deprecated('3.1', alternative=':meth:`Header.totextfile`')
    def toTxtFile(self, fileobj, clobber=False):
        """
        Output the header parameters to a file in ASCII format.

        Parameters
        ----------
        fileobj : file path, file object or file-like object
            Output header parameters file.

        clobber : bool
            When `True`, overwrite the output file if it exists.
        """

        self.tofile(fileobj, sep='\n', endcard=False, padding=False,
                    clobber=clobber)

    @deprecated('3.1',
                message='This is equivalent to '
                        '``self.extend(Header.fromtextfile(fileobj), '
                        'update=True, update_first=True)``.  Note that there '
                        'there is no direct equivalent to the '
                        '``replace=True`` option since '
                        ':meth:`Header.fromtextfile` returns a new '
                        ':class:`Header` instance.')
    def fromTxtFile(self, fileobj, replace=False):
        """
        Input the header parameters from an ASCII file.

        The input header cards will be used to update the current
        header.  Therefore, when an input card key matches a card key
        that already exists in the header, that card will be updated
        in place.  Any input cards that do not already exist in the
        header will be added.  Cards will not be deleted from the
        header.

        Parameters
        ----------
        fileobj : file path, file object or file-like object
            Input header parameters file.

        replace : bool, optional
            When `True`, indicates that the entire header should be
            replaced with the contents of the ASCII file instead of
            just updating the current header.
        """

        input_header = Header.fromfile(fileobj, sep='\n', endcard=False,
                                       padding=False)

        if replace:
            self.clear()
        prev_key = 0

        for card in input_header.cards:
            card.verify('silentfix')

            if card.keyword == 'SIMPLE':
                if self.get('XTENSION'):
                    del self.ascard['XTENSION']

                self.set(card.keyword, card.value, card.comment, before=0)
                prev_key = 0
            elif card.keyword == 'XTENSION':
                if self.get('SIMPLE'):
                    del self.ascard['SIMPLE']

                self.set(card.keyword, card.value, card.comment, before=0)
                prev_key = 0
            elif card.keyword in Card._commentary_keywords:
                if (not replace and
                        not (card.keyword == '' and card.value == '')):
                    # Don't add duplicate commentary cards (though completely
                    # blank cards are allowed to be duplicated)
                    for idx, c in enumerate(self.cards):
                        if c.keyword == card.keyword and c.value == card.value:
                            break
                    else:
                        self.set(card.keyword, card.value, after=prev_key)
                        prev_key += 1
                else:
                    self.set(card.keyword, card.value, after=prev_key)
                    prev_key += 1
            else:
                self.set(card.keyword, card.value, card.comment,
                         after=prev_key)
                prev_key += 1


class _CardAccessor(object):
    """
    This is a generic class for wrapping a Header in such a way that you can
    use the header's slice/filtering capabilities to return a subset of cards
    and do something of them.

    This is sort of the opposite notion of the old CardList class--whereas
    Header used to use CardList to get lists of cards, this uses Header to get
    lists of cards.
    """

    # TODO: Consider giving this dict/list methods like Header itself
    def __init__(self, header):
        self._header = header

    def __repr__(self):
        return '\n'.join(repr(c) for c in self._header._cards)

    def __len__(self):
        return len([c for c in self])

    def __eq__(self, other):
        if isiterable(other):
            for a, b in itertools.izip(self, other):
                if a != b:
                    return False
            else:
                return True
        return False

    def __getitem__(self, item):
        if isinstance(item, slice) or self._header._haswildcard(item):
            return self.__class__(self._header[item])

        idx = self._header._cardindex(item)
        return self._header._cards[idx]

    def _setslice(self, item, value):
        """
        Helper for implementing __setitem__ on _CardAccessor subclasses; slices
        should always be handled in this same way.
        """

        if isinstance(item, slice) or self._header._haswildcard(item):
            if isinstance(item, slice):
                indices = xrange(*item.indices(len(self)))
            else:
                indices = self._header._wildcardmatch(item)
            if isinstance(value, basestring) or not isiterable(value):
                value = itertools.repeat(value, len(indices))
            for idx, val in itertools.izip(indices, value):
                self[idx] = val
            return True
        return False


class _HeaderComments(_CardAccessor):
    """
    A class used internally by the Header class for the Header.comments
    attribute access.

    This object can be used to display all the keyword comments in the Header,
    or look up the comments on specific keywords.  It allows all the same forms
    of keyword lookup as the Header class itself, but returns comments instead
    of values.
    """

    def __repr__(self):
        """Returns a simple list of all keywords and their comments."""

        keyword_length = KEYWORD_LENGTH
        for card in self._header._cards:
            keyword_length = max(keyword_length, len(card.keyword))
        return '\n'.join('%*s  %s' % (keyword_length, c.keyword, c.comment)
                         for c in self._header._cards)

    def __getitem__(self, item):
        """
        Slices and filter strings return a new _HeaderComments containing the
        returned cards.  Otherwise the comment of a single card is returned.
        """

        item = super(_HeaderComments, self).__getitem__(item)
        if isinstance(item, _HeaderComments):
            return item
        return item.comment

    def __setitem__(self, item, comment):
        """
        Set/update the comment on specified card or cards.

        Slice/filter updates work similarly to how Header.__setitem__ works.
        """

        if self._setslice(item, comment):
            return

        # In this case, key/index errors should be raised; don't update
        # comments of nonexistent cards
        idx = self._header._cardindex(item)
        value = self._header[idx]
        self._header[idx] = (value, comment)


class _HeaderCommentaryCards(_CardAccessor):
    def __init__(self, header, keyword=''):
        super(_HeaderCommentaryCards, self).__init__(header)
        self._keyword = keyword
        self._count = self._header.count(self._keyword)
        self._indices = slice(self._count).indices(self._count)

    def __repr__(self):
        return '\n'.join(self)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            n = self.__class__(self._header, self._keyword)
            n._indices = idx.indices(self._count)
            return n
        elif not isinstance(idx, int):
            raise ValueError('%s index must be an integer' % self._keyword)

        idx = range(*self._indices)[idx]
        return self._header[(self._keyword, idx)]

    def __setitem__(self, item, value):
        """
        Set the value of a specified commentary card or cards.

        Slice/filter updates work similarly to how Header.__setitem__ works.
        """

        if self._setslice(item, value):
            return

        # In this case, key/index errors should be raised; don't update
        # comments of nonexistent cards
        self._header[(self._keyword, item)] = value


def _is_astropy_internal():
    """
    Returns True if the stack frame this is called from is in code internal to
    the the astropy package.

    This is used in a few places where hacks are employed for backwards
    compatibility with the old header API, but where we want to avoid using
    those hacks internally.
    """

    calling_mod = inspect.getmodule(sys._getframe(2))
    return calling_mod and calling_mod.__name__.startswith('astropy.')


def _block_size(sep):
    """
    Determine the size of a FITS header block if a non-blank separator is used
    between cards.
    """

    return BLOCK_SIZE + (len(sep) * (BLOCK_SIZE // Card.length - 1))


def _pad_length(stringlen):
    """Bytes needed to pad the input stringlen to the next FITS block."""

    return (BLOCK_SIZE - (stringlen % BLOCK_SIZE)) % BLOCK_SIZE
