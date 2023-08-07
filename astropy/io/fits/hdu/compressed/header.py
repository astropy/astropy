# Licensed under a 3-clause BSD style license - see PYFITS.rst

import re
import warnings
from contextlib import suppress

from astropy.io.fits.card import Card
from astropy.io.fits.column import KEYWORD_NAMES as TABLE_KEYWORD_NAMES
from astropy.io.fits.column import TDEF_RE
from astropy.io.fits.header import Header

__all__ = ['CompImageHeader']


class CompImageHeader(Header):
    """
    Header object for compressed image HDUs designed to keep the compression
    header and the underlying image header properly synchronized.

    This essentially wraps the image header, so that all values are read from
    and written to the image header.  However, updates to the image header will
    also update the table header where appropriate.

    Note that if no image header is passed in, the code will instantiate a
    regular `~astropy.io.fits.Header`.
    """

    # TODO: The difficulty of implementing this screams a need to rewrite this
    # module

    _keyword_remaps = {
        "SIMPLE": "ZSIMPLE",
        "XTENSION": "ZTENSION",
        "BITPIX": "ZBITPIX",
        "NAXIS": "ZNAXIS",
        "EXTEND": "ZEXTEND",
        "BLOCKED": "ZBLOCKED",
        "PCOUNT": "ZPCOUNT",
        "GCOUNT": "ZGCOUNT",
        "CHECKSUM": "ZHECKSUM",
        "DATASUM": "ZDATASUM",
    }

    _zdef_re = re.compile(r"(?P<label>^[Zz][a-zA-Z]*)(?P<num>[1-9][0-9 ]*$)?")
    _compression_keywords = set(_keyword_remaps.values()).union(
        ["ZIMAGE", "ZCMPTYPE", "ZMASKCMP", "ZQUANTIZ", "ZDITHER0"]
    )
    _indexed_compression_keywords = {"ZNAXIS", "ZTILE", "ZNAME", "ZVAL"}
    # TODO: Once it place it should be possible to manage some of this through
    # the schema system, but it's not quite ready for that yet.  Also it still
    # makes more sense to change CompImageHDU to subclass ImageHDU :/

    def __new__(cls, table_header, image_header=None):
        # 2019-09-14 (MHvK): No point wrapping anything if no image_header is
        # given.  This happens if __getitem__ and copy are called - our super
        # class will aim to initialize a new, possibly partially filled
        # header, but we cannot usefully deal with that.
        # TODO: the above suggests strongly we should *not* subclass from
        # Header.  See also comment above about the need for reorganization.
        if image_header is None:
            return Header(table_header)
        else:
            return super().__new__(cls)

    def __init__(self, table_header, image_header):
        self._cards = image_header._cards
        self._keyword_indices = image_header._keyword_indices
        self._rvkc_indices = image_header._rvkc_indices
        self._modified = image_header._modified
        self._table_header = table_header

    # We need to override and Header methods that can modify the header, and
    # ensure that they sync with the underlying _table_header

    def __setitem__(self, key, value):
        # This isn't pretty, but if the `key` is either an int or a tuple we
        # need to figure out what keyword name that maps to before doing
        # anything else; these checks will be repeated later in the
        # super().__setitem__ call but I don't see another way around it
        # without some major refactoring
        if self._set_slice(key, value, self):
            return

        if isinstance(key, int):
            keyword, index = self._keyword_from_index(key)
        elif isinstance(key, tuple):
            keyword, index = key
        else:
            # We don't want to specify and index otherwise, because that will
            # break the behavior for new keywords and for commentary keywords
            keyword, index = key, None

        if self._is_reserved_keyword(keyword):
            return

        super().__setitem__(key, value)

        if index is not None:
            remapped_keyword = self._remap_keyword(keyword)
            self._table_header[remapped_keyword, index] = value
        # Else this will pass through to ._update

    def __delitem__(self, key):
        if isinstance(key, slice) or self._haswildcard(key):
            # If given a slice pass that on to the superclass and bail out
            # early; we only want to make updates to _table_header when given
            # a key specifying a single keyword
            return super().__delitem__(key)

        if isinstance(key, int):
            keyword, index = self._keyword_from_index(key)
        elif isinstance(key, tuple):
            keyword, index = key
        else:
            keyword, index = key, None

        if key not in self:
            raise KeyError(f"Keyword {key!r} not found.")

        super().__delitem__(key)

        remapped_keyword = self._remap_keyword(keyword)

        if remapped_keyword in self._table_header:
            if index is not None:
                del self._table_header[(remapped_keyword, index)]
            else:
                del self._table_header[remapped_keyword]

    def append(self, card=None, useblanks=True, bottom=False, end=False):
        # This logic unfortunately needs to be duplicated from the base class
        # in order to determine the keyword
        if isinstance(card, str):
            card = Card(card)
        elif isinstance(card, tuple):
            card = Card(*card)
        elif card is None:
            card = Card()
        elif not isinstance(card, Card):
            raise ValueError(
                "The value appended to a Header must be either a keyword or "
                "(keyword, value, [comment]) tuple; got: {!r}".format(card)
            )

        if self._is_reserved_keyword(card.keyword):
            return

        super().append(card=card, useblanks=useblanks, bottom=bottom, end=end)

        remapped_keyword = self._remap_keyword(card.keyword)

        # card.keyword strips the HIERARCH if present so this must be added
        # back to avoid a warning.
        if str(card).startswith("HIERARCH ") and not remapped_keyword.startswith(
            "HIERARCH "
        ):
            remapped_keyword = "HIERARCH " + remapped_keyword

        card = Card(remapped_keyword, card.value, card.comment)

        # Here we disable the use of blank cards, because the call above to
        # Header.append may have already deleted a blank card in the table
        # header, thanks to inheritance: Header.append calls 'del self[-1]'
        # to delete a blank card, which calls CompImageHeader.__deltitem__,
        # which deletes the blank card both in the image and the table headers!
        self._table_header.append(card=card, useblanks=False, bottom=bottom, end=end)

    def insert(self, key, card, useblanks=True, after=False):
        if isinstance(key, int):
            # Determine condition to pass through to append
            if after:
                if key == -1:
                    key = len(self._cards)
                else:
                    key += 1

            if key >= len(self._cards):
                self.append(card, end=True)
                return

        if isinstance(card, str):
            card = Card(card)
        elif isinstance(card, tuple):
            card = Card(*card)
        elif not isinstance(card, Card):
            raise ValueError(
                "The value inserted into a Header must be either a keyword or "
                "(keyword, value, [comment]) tuple; got: {!r}".format(card)
            )

        if self._is_reserved_keyword(card.keyword):
            return

        # Now the tricky part is to determine where to insert in the table
        # header.  If given a numerical index we need to map that to the
        # corresponding index in the table header.  Although rare, there may be
        # cases where there is no mapping in which case we just try the same
        # index
        # NOTE: It is crucial that remapped_index in particular is figured out
        # before the image header is modified
        remapped_index = self._remap_index(key)
        remapped_keyword = self._remap_keyword(card.keyword)

        super().insert(key, card, useblanks=useblanks, after=after)

        card = Card(remapped_keyword, card.value, card.comment)

        # Here we disable the use of blank cards, because the call above to
        # Header.insert may have already deleted a blank card in the table
        # header, thanks to inheritance: Header.insert calls 'del self[-1]'
        # to delete a blank card, which calls CompImageHeader.__delitem__,
        # which deletes the blank card both in the image and the table headers!
        self._table_header.insert(remapped_index, card, useblanks=False, after=after)

    def _update(self, card):
        keyword = card[0]

        if self._is_reserved_keyword(keyword):
            return

        super()._update(card)

        if keyword in Card._commentary_keywords:
            # Otherwise this will result in a duplicate insertion
            return

        remapped_keyword = self._remap_keyword(keyword)
        self._table_header._update((remapped_keyword,) + card[1:])

    # Last piece needed (I think) for synchronizing with the real header
    # This one is tricky since _relativeinsert calls insert
    def _relativeinsert(self, card, before=None, after=None, replace=False):
        keyword = card[0]

        if self._is_reserved_keyword(keyword):
            return

        # Now we have to figure out how to remap 'before' and 'after'
        if before is None:
            if isinstance(after, int):
                remapped_after = self._remap_index(after)
            else:
                remapped_after = self._remap_keyword(after)
            remapped_before = None
        else:
            if isinstance(before, int):
                remapped_before = self._remap_index(before)
            else:
                remapped_before = self._remap_keyword(before)
            remapped_after = None

        super()._relativeinsert(card, before=before, after=after, replace=replace)

        remapped_keyword = self._remap_keyword(keyword)

        card = Card(remapped_keyword, card[1], card[2])
        self._table_header._relativeinsert(
            card, before=remapped_before, after=remapped_after, replace=replace
        )

    @classmethod
    def _is_reserved_keyword(cls, keyword, warn=True):
        msg = (
            "Keyword {!r} is reserved for use by the FITS Tiled Image "
            "Convention and will not be stored in the header for the "
            "image being compressed.".format(keyword)
        )

        if keyword == "TFIELDS":
            if warn:
                warnings.warn(msg)
            return True

        m = TDEF_RE.match(keyword)

        if m and m.group("label").upper() in TABLE_KEYWORD_NAMES:
            if warn:
                warnings.warn(msg)
            return True

        m = cls._zdef_re.match(keyword)

        if m:
            label = m.group("label").upper()
            num = m.group("num")
            if num is not None and label in cls._indexed_compression_keywords:
                if warn:
                    warnings.warn(msg)
                return True
            elif label in cls._compression_keywords:
                if warn:
                    warnings.warn(msg)
                return True

        return False

    @classmethod
    def _remap_keyword(cls, keyword):
        # Given a keyword that one might set on an image, remap that keyword to
        # the name used for it in the COMPRESSED HDU header
        # This is mostly just a lookup in _keyword_remaps, but needs handling
        # for NAXISn keywords

        is_naxisn = False
        if keyword[:5] == "NAXIS":
            with suppress(ValueError):
                index = int(keyword[5:])
                is_naxisn = index > 0

        if is_naxisn:
            return f"ZNAXIS{index}"

        # If the keyword does not need to be remapped then just return the
        # original keyword
        return cls._keyword_remaps.get(keyword, keyword)

    def _remap_index(self, idx):
        # Given an integer index into this header, map that to the index in the
        # table header for the same card.  If the card doesn't exist in the
        # table header (generally should *not* be the case) this will just
        # return the same index
        # This *does* also accept a keyword or (keyword, repeat) tuple and
        # obtains the associated numerical index with self._cardindex
        if not isinstance(idx, int):
            idx = self._cardindex(idx)

        keyword, repeat = self._keyword_from_index(idx)
        remapped_insert_keyword = self._remap_keyword(keyword)

        with suppress(IndexError, KeyError):
            idx = self._table_header._cardindex((remapped_insert_keyword, repeat))

        return idx

    def clear(self):
        """
        Remove all cards from the header.
        """
        self._table_header.clear()
        super().clear()
