# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools
import re
import warnings
from contextlib import suppress

from astropy.io.fits.card import Card
from astropy.io.fits.column import KEYWORD_NAMES as TABLE_KEYWORD_NAMES
from astropy.io.fits.column import TDEF_RE, ColDefs, Column
from astropy.io.fits.header import Header
from astropy.utils.exceptions import AstropyUserWarning

from .settings import (
    CMTYPE_ALIASES,
    COMPRESSION_TYPES,
    DEFAULT_BLOCK_SIZE,
    DEFAULT_BYTE_PIX,
    DEFAULT_COMPRESSION_TYPE,
    DEFAULT_DITHER_SEED,
    DEFAULT_HCOMP_SCALE,
    DEFAULT_HCOMP_SMOOTH,
    DEFAULT_QUANTIZE_LEVEL,
    DEFAULT_QUANTIZE_METHOD,
    NO_DITHER,
    QUANTIZE_METHOD_NAMES,
    SUBTRACTIVE_DITHER_1,
    SUBTRACTIVE_DITHER_2,
)
from .utils import _validate_tile_shape

__all__ = [
    "CompImageHeader",
    "_bintable_header_to_image_header",
    "_image_header_to_bintable_header_and_coldefs",
]


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
                f"(keyword, value, [comment]) tuple; got: {card!r}"
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
                f"(keyword, value, [comment]) tuple; got: {card!r}"
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
            f"Keyword {keyword!r} is reserved for use by the FITS Tiled Image "
            "Convention and will not be stored in the header for the "
            "image being compressed."
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


def _bintable_header_to_image_header(bintable_header):
    # Start with a copy of the table header.
    image_header = bintable_header.copy()

    # Delete cards that are related to the table.  And move
    # the values of those cards that relate to the image from
    # their corresponding table cards.  These include
    # ZBITPIX -> BITPIX, ZNAXIS -> NAXIS, and ZNAXISn -> NAXISn.
    # (Note: Used set here instead of list in case there are any duplicate
    # keywords, which there may be in some pathological cases:
    # https://github.com/astropy/astropy/issues/2750
    for keyword in set(image_header):
        if CompImageHeader._is_reserved_keyword(keyword, warn=False):
            del image_header[keyword]

    hcomments = bintable_header.comments

    if "ZSIMPLE" in bintable_header:
        image_header.set(
            "SIMPLE", bintable_header["ZSIMPLE"], hcomments["ZSIMPLE"], before=0
        )
        del image_header["XTENSION"]
    elif "ZTENSION" in bintable_header:
        if bintable_header["ZTENSION"] != "IMAGE":
            warnings.warn(
                "ZTENSION keyword in compressed extension != 'IMAGE'",
                AstropyUserWarning,
            )
        image_header.set("XTENSION", "IMAGE", hcomments["ZTENSION"], before=0)
    else:
        image_header.set("XTENSION", "IMAGE", before=0)

    image_header.set(
        "BITPIX", bintable_header["ZBITPIX"], hcomments["ZBITPIX"], before=1
    )

    image_header.set("NAXIS", bintable_header["ZNAXIS"], hcomments["ZNAXIS"], before=2)

    last_naxis = "NAXIS"
    for idx in range(image_header["NAXIS"]):
        znaxis = "ZNAXIS" + str(idx + 1)
        naxis = znaxis[1:]
        image_header.set(
            naxis, bintable_header[znaxis], hcomments[znaxis], after=last_naxis
        )
        last_naxis = naxis

    # Delete any other spurious NAXISn keywords:
    naxis = image_header["NAXIS"]
    for keyword in list(image_header["NAXIS?*"]):
        try:
            n = int(keyword[5:])
        except Exception:
            continue

        if n > naxis:
            del image_header[keyword]

    # Although PCOUNT and GCOUNT are considered mandatory for IMAGE HDUs,
    # ZPCOUNT and ZGCOUNT are optional, probably because for IMAGE HDUs
    # their values are always 0 and 1 respectively
    if "ZPCOUNT" in bintable_header:
        image_header.set(
            "PCOUNT",
            bintable_header["ZPCOUNT"],
            hcomments["ZPCOUNT"],
            after=last_naxis,
        )
    else:
        image_header.set("PCOUNT", 0, after=last_naxis)

    if "ZGCOUNT" in bintable_header:
        image_header.set(
            "GCOUNT", bintable_header["ZGCOUNT"], hcomments["ZGCOUNT"], after="PCOUNT"
        )
    else:
        image_header.set("GCOUNT", 1, after="PCOUNT")

    if "ZEXTEND" in bintable_header:
        image_header.set("EXTEND", bintable_header["ZEXTEND"], hcomments["ZEXTEND"])

    if "ZBLOCKED" in bintable_header:
        image_header.set("BLOCKED", bintable_header["ZBLOCKED"], hcomments["ZBLOCKED"])

    # Move the ZHECKSUM and ZDATASUM cards to the image header
    # as CHECKSUM and DATASUM
    if "ZHECKSUM" in bintable_header:
        image_header.set("CHECKSUM", bintable_header["ZHECKSUM"], hcomments["ZHECKSUM"])

    if "ZDATASUM" in bintable_header:
        image_header.set("DATASUM", bintable_header["ZDATASUM"], hcomments["ZDATASUM"])

    # Remove the EXTNAME card if the value in the table header
    # is the default value of COMPRESSED_IMAGE.
    if "EXTNAME" in image_header and image_header["EXTNAME"] == "COMPRESSED_IMAGE":
        del image_header["EXTNAME"]

    # Remove the PCOUNT GCOUNT cards if the uncompressed header is
    # from a primary HDU
    if "SIMPLE" in image_header:
        del image_header["PCOUNT"]
        del image_header["GCOUNT"]

    # Look to see if there are any blank cards in the table
    # header.  If there are, there should be the same number
    # of blank cards in the image header.  Add blank cards to
    # the image header to make it so.
    table_blanks = bintable_header._countblanks()
    image_blanks = image_header._countblanks()

    for _ in range(table_blanks - image_blanks):
        image_header.append()

    # Create the CompImageHeader that syncs with the table header
    return CompImageHeader(bintable_header, image_header)


def _image_header_to_bintable_header_and_coldefs(
    image_header,
    generated_image_header,
    bintable_header,
    name=None,
    huge_hdu=False,
    compression_type=None,
    tile_shape=None,
    hcomp_scale=None,
    hcomp_smooth=None,
    quantize_level=None,
    quantize_method=None,
    dither_seed=None,
    axes=None,
    generate_dither_seed=None,
):
    # NOTE: image_header is the header that a user would see as the image
    # header which they might have set things like BSCALE and BZERO on, or
    # added history or comments to. Whereas generated_image_header is the
    # image header as converted/generated from the existing binary table HDU.

    # Update the extension name in the table header
    if not name and "EXTNAME" not in bintable_header:
        # Do not sync this with the image header since the default
        # name is specific to the table header.
        bintable_header.set(
            "EXTNAME",
            "COMPRESSED_IMAGE",
            "name of this binary table extension",
            after="TFIELDS",
        )

    # Set the compression type in the table header.
    if compression_type:
        if compression_type not in COMPRESSION_TYPES:
            warnings.warn(
                "Unknown compression type provided (supported are {}). "
                "Default ({}) compression will be used.".format(
                    ", ".join(map(repr, COMPRESSION_TYPES)),
                    DEFAULT_COMPRESSION_TYPE,
                ),
                AstropyUserWarning,
            )
            compression_type = DEFAULT_COMPRESSION_TYPE

        bintable_header.set(
            "ZCMPTYPE", compression_type, "compression algorithm", after="TFIELDS"
        )
    else:
        compression_type = CMTYPE_ALIASES.get(compression_type, compression_type)

    # If the input image header had BSCALE/BZERO cards, then insert
    # them in the table header.

    if image_header:
        bzero = image_header.get("BZERO", 0.0)
        bscale = image_header.get("BSCALE", 1.0)
        after_keyword = "EXTNAME"

        if bscale != 1.0:
            bintable_header.set("BSCALE", bscale, after=after_keyword)
            after_keyword = "BSCALE"

        if bzero != 0.0:
            bintable_header.set("BZERO", bzero, after=after_keyword)

    try:
        bitpix_comment = image_header.comments["BITPIX"]
    except (AttributeError, KeyError):
        bitpix_comment = "data type of original image"

    try:
        naxis_comment = image_header.comments["NAXIS"]
    except (AttributeError, KeyError):
        naxis_comment = "dimension of original image"

    # Set the label for the first column in the table

    bintable_header.set(
        "TTYPE1", "COMPRESSED_DATA", "label for field 1", after="TFIELDS"
    )

    # Set the data format for the first column.  It is dependent
    # on the requested compression type.

    if compression_type == "PLIO_1":
        tform1 = "1QI" if huge_hdu else "1PI"
    else:
        tform1 = "1QB" if huge_hdu else "1PB"

    bintable_header.set(
        "TFORM1",
        tform1,
        "data format of field: variable length array",
        after="TTYPE1",
    )

    # Create the first column for the table.  This column holds the
    # compressed data.
    col1 = Column(name=bintable_header["TTYPE1"], format=tform1)

    # Create the additional columns required for floating point
    # data and calculate the width of the output table.

    zbitpix = generated_image_header["BITPIX"]

    if zbitpix < 0 and quantize_level != 0.0:
        # floating point image has 'COMPRESSED_DATA',
        # 'UNCOMPRESSED_DATA', 'ZSCALE', and 'ZZERO' columns (unless using
        # lossless compression, per CFITSIO)
        ncols = 4

        # CFITSIO 3.28 and up automatically use the GZIP_COMPRESSED_DATA
        # store floating point data that couldn't be quantized, instead
        # of the UNCOMPRESSED_DATA column.  There's no way to control
        # this behavior so the only way to determine which behavior will
        # be employed is via the CFITSIO version

        ttype2 = "GZIP_COMPRESSED_DATA"
        # The required format for the GZIP_COMPRESSED_DATA is actually
        # missing from the standard docs, but CFITSIO suggests it
        # should be 1PB, which is logical.
        tform2 = "1QB" if huge_hdu else "1PB"

        # Set up the second column for the table that will hold any
        # uncompressable data.
        bintable_header.set("TTYPE2", ttype2, "label for field 2", after="TFORM1")

        bintable_header.set(
            "TFORM2",
            tform2,
            "data format of field: variable length array",
            after="TTYPE2",
        )

        col2 = Column(name=ttype2, format=tform2)

        # Set up the third column for the table that will hold
        # the scale values for quantized data.
        bintable_header.set("TTYPE3", "ZSCALE", "label for field 3", after="TFORM2")
        bintable_header.set(
            "TFORM3", "1D", "data format of field: 8-byte DOUBLE", after="TTYPE3"
        )
        col3 = Column(name=bintable_header["TTYPE3"], format=bintable_header["TFORM3"])

        # Set up the fourth column for the table that will hold
        # the zero values for the quantized data.
        bintable_header.set("TTYPE4", "ZZERO", "label for field 4", after="TFORM3")
        bintable_header.set(
            "TFORM4", "1D", "data format of field: 8-byte DOUBLE", after="TTYPE4"
        )
        after = "TFORM4"
        col4 = Column(name=bintable_header["TTYPE4"], format=bintable_header["TFORM4"])

        # Create the ColDefs object for the table
        cols = ColDefs([col1, col2, col3, col4])
    else:
        # default table has just one 'COMPRESSED_DATA' column
        ncols = 1
        after = "TFORM1"

        # remove any header cards for the additional columns that
        # may be left over from the previous data
        to_remove = ["TTYPE2", "TFORM2", "TTYPE3", "TFORM3", "TTYPE4", "TFORM4"]

        for k in to_remove:
            try:
                del bintable_header[k]
            except KeyError:
                pass

        # Create the ColDefs object for the table
        cols = ColDefs([col1])

    # Update the table header with the width of the table, the
    # number of fields in the table, the indicator for a compressed
    # image HDU, the data type of the image data and the number of
    # dimensions in the image data array.
    bintable_header.set("NAXIS1", cols.dtype.itemsize, "width of table in bytes")
    bintable_header.set(
        "TFIELDS", ncols, "number of fields in each row", after="GCOUNT"
    )
    bintable_header.set(
        "ZIMAGE", True, "extension contains compressed image", after=after
    )
    bintable_header.set("ZBITPIX", zbitpix, bitpix_comment, after="ZIMAGE")
    bintable_header.set(
        "ZNAXIS", generated_image_header["NAXIS"], naxis_comment, after="ZBITPIX"
    )

    # Strip the table header of all the ZNAZISn and ZTILEn keywords
    # that may be left over from the previous data

    for idx in itertools.count(1):
        try:
            del bintable_header["ZNAXIS" + str(idx)]
            del bintable_header["ZTILE" + str(idx)]
        except KeyError:
            break

    # Verify that any input tile size parameter is the appropriate
    # size to match the HDU's data.

    tile_shape = _validate_tile_shape(
        tile_shape=tile_shape,
        compression_type=compression_type,
        image_header=generated_image_header,
    )

    # Set up locations for writing the next cards in the header.
    last_znaxis = "ZNAXIS"

    if generated_image_header["NAXIS"] > 0:
        after1 = "ZNAXIS1"
    else:
        after1 = "ZNAXIS"

    # Calculate the number of rows in the output table and
    # write the ZNAXISn and ZTILEn cards to the table header.
    nrows = 0

    for idx, axis in enumerate(axes):
        naxis = "NAXIS" + str(idx + 1)
        znaxis = "ZNAXIS" + str(idx + 1)
        ztile = "ZTILE" + str(idx + 1)

        ts = tile_shape[len(axes) - 1 - idx]

        if not nrows:
            nrows = (axis - 1) // ts + 1
        else:
            nrows *= (axis - 1) // ts + 1

        if image_header and naxis in image_header:
            bintable_header.set(
                znaxis, axis, image_header.comments[naxis], after=last_znaxis
            )
        else:
            bintable_header.set(
                znaxis, axis, "length of original image axis", after=last_znaxis
            )

        bintable_header.set(ztile, ts, "size of tiles to be compressed", after=after1)
        last_znaxis = znaxis
        after1 = ztile

    # Set the NAXIS2 header card in the table hdu to the number of
    # rows in the table.
    bintable_header.set("NAXIS2", nrows, "number of rows in table")

    # Set the compression parameters in the table header.

    # First, setup the values to be used for the compression parameters
    # in case none were passed in.  This will be either the value
    # already in the table header for that parameter or the default
    # value.
    for idx in itertools.count(1):
        zname = "ZNAME" + str(idx)
        if zname not in bintable_header:
            break
        zval = "ZVAL" + str(idx)
        if bintable_header[zname] == "NOISEBIT":
            if quantize_level is None:
                quantize_level = bintable_header[zval]
        if bintable_header[zname] == "SCALE   ":
            if hcomp_scale is None:
                hcomp_scale = bintable_header[zval]
        if bintable_header[zname] == "SMOOTH  ":
            if hcomp_smooth is None:
                hcomp_smooth = bintable_header[zval]

    if quantize_level is None:
        quantize_level = DEFAULT_QUANTIZE_LEVEL

    if hcomp_scale is None:
        hcomp_scale = DEFAULT_HCOMP_SCALE

    if hcomp_smooth is None:
        hcomp_smooth = DEFAULT_HCOMP_SMOOTH

    # Next, strip the table header of all the ZNAMEn and ZVALn keywords
    # that may be left over from the previous data
    for idx in itertools.count(1):
        zname = "ZNAME" + str(idx)
        if zname not in bintable_header:
            break
        zval = "ZVAL" + str(idx)
        del bintable_header[zname]
        del bintable_header[zval]

    # Finally, put the appropriate keywords back based on the
    # compression type.

    after_keyword = "ZCMPTYPE"
    idx = 1

    if compression_type == "RICE_1":
        bintable_header.set(
            "ZNAME1", "BLOCKSIZE", "compression block size", after=after_keyword
        )
        bintable_header.set(
            "ZVAL1", DEFAULT_BLOCK_SIZE, "pixels per block", after="ZNAME1"
        )

        bintable_header.set(
            "ZNAME2", "BYTEPIX", "bytes per pixel (1, 2, 4, or 8)", after="ZVAL1"
        )

        if bintable_header["ZBITPIX"] == 8:
            bytepix = 1
        elif bintable_header["ZBITPIX"] == 16:
            bytepix = 2
        else:
            bytepix = DEFAULT_BYTE_PIX

        bintable_header.set(
            "ZVAL2", bytepix, "bytes per pixel (1, 2, 4, or 8)", after="ZNAME2"
        )
        after_keyword = "ZVAL2"
        idx = 3
    elif compression_type == "HCOMPRESS_1":
        bintable_header.set(
            "ZNAME1", "SCALE", "HCOMPRESS scale factor", after=after_keyword
        )
        bintable_header.set(
            "ZVAL1", hcomp_scale, "HCOMPRESS scale factor", after="ZNAME1"
        )
        bintable_header.set(
            "ZNAME2", "SMOOTH", "HCOMPRESS smooth option", after="ZVAL1"
        )
        bintable_header.set(
            "ZVAL2", hcomp_smooth, "HCOMPRESS smooth option", after="ZNAME2"
        )
        after_keyword = "ZVAL2"
        idx = 3

    if generated_image_header["BITPIX"] < 0:  # floating point image
        bintable_header.set(
            "ZNAME" + str(idx),
            "NOISEBIT",
            "floating point quantization level",
            after=after_keyword,
        )
        bintable_header.set(
            "ZVAL" + str(idx),
            quantize_level,
            "floating point quantization level",
            after="ZNAME" + str(idx),
        )

        # Add the dither method and seed
        if quantize_method:
            if quantize_method not in [
                NO_DITHER,
                SUBTRACTIVE_DITHER_1,
                SUBTRACTIVE_DITHER_2,
            ]:
                name = QUANTIZE_METHOD_NAMES[DEFAULT_QUANTIZE_METHOD]
                warnings.warn(
                    "Unknown quantization method provided.  "
                    f"Default method ({name}) used."
                )
                quantize_method = DEFAULT_QUANTIZE_METHOD

            if quantize_method == NO_DITHER:
                zquantiz_comment = "No dithering during quantization"
            else:
                zquantiz_comment = "Pixel Quantization Algorithm"

            bintable_header.set(
                "ZQUANTIZ",
                QUANTIZE_METHOD_NAMES[quantize_method],
                zquantiz_comment,
                after="ZVAL" + str(idx),
            )
        else:
            # If the ZQUANTIZ keyword is missing the default is to assume
            # no dithering, rather than whatever DEFAULT_QUANTIZE_METHOD
            # is set to
            quantize_method = bintable_header.get("ZQUANTIZ", NO_DITHER)

            if isinstance(quantize_method, str):
                for k, v in QUANTIZE_METHOD_NAMES.items():
                    if v.upper() == quantize_method:
                        quantize_method = k
                        break
                else:
                    quantize_method = NO_DITHER

        if quantize_method == NO_DITHER:
            if "ZDITHER0" in bintable_header:
                # If dithering isn't being used then there's no reason to
                # keep the ZDITHER0 keyword
                del bintable_header["ZDITHER0"]
        else:
            if dither_seed:
                dither_seed = generate_dither_seed(dither_seed)
            elif "ZDITHER0" in bintable_header:
                dither_seed = bintable_header["ZDITHER0"]
            else:
                dither_seed = generate_dither_seed(DEFAULT_DITHER_SEED)

            bintable_header.set(
                "ZDITHER0",
                dither_seed,
                "dithering offset when quantizing floats",
                after="ZQUANTIZ",
            )

    if image_header:
        # Move SIMPLE card from the image header to the
        # table header as ZSIMPLE card.

        if "SIMPLE" in image_header:
            bintable_header.set(
                "ZSIMPLE",
                image_header["SIMPLE"],
                image_header.comments["SIMPLE"],
                before="ZBITPIX",
            )

        # Move EXTEND card from the image header to the
        # table header as ZEXTEND card.

        if "EXTEND" in image_header:
            bintable_header.set(
                "ZEXTEND", image_header["EXTEND"], image_header.comments["EXTEND"]
            )

        # Move BLOCKED card from the image header to the
        # table header as ZBLOCKED card.

        if "BLOCKED" in image_header:
            bintable_header.set(
                "ZBLOCKED",
                image_header["BLOCKED"],
                image_header.comments["BLOCKED"],
            )

        # Move XTENSION card from the image header to the
        # table header as ZTENSION card.

        # Since we only handle compressed IMAGEs, ZTENSION should
        # always be IMAGE, even if the caller has passed in a header
        # for some other type of extension.
        if "XTENSION" in image_header:
            bintable_header.set(
                "ZTENSION",
                "IMAGE",
                image_header.comments["XTENSION"],
                before="ZBITPIX",
            )

        # Move PCOUNT and GCOUNT cards from image header to the table
        # header as ZPCOUNT and ZGCOUNT cards.

        if "PCOUNT" in image_header:
            bintable_header.set(
                "ZPCOUNT",
                image_header["PCOUNT"],
                image_header.comments["PCOUNT"],
                after=last_znaxis,
            )

        if "GCOUNT" in image_header:
            bintable_header.set(
                "ZGCOUNT",
                image_header["GCOUNT"],
                image_header.comments["GCOUNT"],
                after="ZPCOUNT",
            )

        # Move CHECKSUM and DATASUM cards from the image header to the
        # table header as XHECKSUM and XDATASUM cards.

        if "CHECKSUM" in image_header:
            bintable_header.set(
                "ZHECKSUM",
                image_header["CHECKSUM"],
                image_header.comments["CHECKSUM"],
            )

        if "DATASUM" in image_header:
            bintable_header.set(
                "ZDATASUM",
                image_header["DATASUM"],
                image_header.comments["DATASUM"],
            )
    else:
        # Move XTENSION card from the image header to the
        # table header as ZTENSION card.

        # Since we only handle compressed IMAGEs, ZTENSION should
        # always be IMAGE, even if the caller has passed in a header
        # for some other type of extension.
        if "XTENSION" in generated_image_header:
            bintable_header.set(
                "ZTENSION",
                "IMAGE",
                generated_image_header.comments["XTENSION"],
                before="ZBITPIX",
            )

        # Move PCOUNT and GCOUNT cards from image header to the table
        # header as ZPCOUNT and ZGCOUNT cards.

        if "PCOUNT" in generated_image_header:
            bintable_header.set(
                "ZPCOUNT",
                generated_image_header["PCOUNT"],
                generated_image_header.comments["PCOUNT"],
                after=last_znaxis,
            )

        if "GCOUNT" in generated_image_header:
            bintable_header.set(
                "ZGCOUNT",
                generated_image_header["GCOUNT"],
                generated_image_header.comments["GCOUNT"],
                after="ZPCOUNT",
            )

    # When we have an image checksum we need to ensure that the same
    # number of blank cards exist in the table header as there were in
    # the image header.  This allows those blank cards to be carried
    # over to the image header when the hdu is uncompressed.

    if "ZHECKSUM" in bintable_header:
        required_blanks = image_header._countblanks()
        image_blanks = generated_image_header._countblanks()
        table_blanks = bintable_header._countblanks()

        for _ in range(required_blanks - image_blanks):
            generated_image_header.append()
            table_blanks += 1

        for _ in range(required_blanks - table_blanks):
            bintable_header.append()

    return bintable_header, cols
