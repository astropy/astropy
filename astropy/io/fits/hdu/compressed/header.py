# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import warnings

from astropy.io.fits.card import Card
from astropy.io.fits.column import KEYWORD_NAMES as TABLE_KEYWORD_NAMES
from astropy.io.fits.column import TDEF_RE, ColDefs, Column
from astropy.io.fits.hdu.compressed.compbintable import _CompBinTableHDU
from astropy.io.fits.header import Header
from astropy.utils.exceptions import AstropyUserWarning

from .settings import (
    CMTYPE_ALIASES,
    COMPRESSION_TYPES,
    DEFAULT_BLOCK_SIZE,
    DEFAULT_BYTE_PIX,
    DEFAULT_COMPRESSION_TYPE,
    DEFAULT_DITHER_SEED,
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
    "_image_header_to_empty_bintable",
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

        super().insert(key, card, useblanks=useblanks, after=after)

    def _update(self, card):
        keyword = card[0]

        if self._is_reserved_keyword(keyword):
            return

        super()._update(card)

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


def _bintable_header_to_image_header(bintable_header):
    # Start with a copy of the table header.
    image_header = bintable_header.copy()

    bscale = image_header.get("BSCALE")
    bzero = image_header.get("BZERO")

    # Strip out special keywords
    image_header.strip()

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

    if bscale:
        image_header['BSCALE'] = bscale
    if bzero:
        image_header['BZERO'] = bzero

    hcomments = bintable_header.comments

    if "ZSIMPLE" in bintable_header:
        image_header.set(
            "SIMPLE", bintable_header["ZSIMPLE"], hcomments["ZSIMPLE"], before=0
        )
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
    return CompImageHeader(image_header)


def _image_header_to_empty_bintable(
    image_header,
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
    bintable = _CompBinTableHDU()

    # NOTE: image_header is the header that a user would see as the image
    # header which they might have set things like BSCALE and BZERO on, or
    # added history or comments to.

    # Update the extension name in the table header
    if name:
        bintable.header["EXTNAME"] = name
    else:
        # Do not sync this with the image header since the default
        # name is specific to the table header.
        bintable.header.set(
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

        bintable.header.set(
            "ZCMPTYPE", compression_type, "compression algorithm", after="TFIELDS"
        )
    else:
        compression_type = CMTYPE_ALIASES.get(compression_type, compression_type)

    # If the input image header had BSCALE/BZERO cards, then insert
    # them in the table header.

    bzero = image_header.get("BZERO", 0.0)
    bscale = image_header.get("BSCALE", 1.0)
    after_keyword = "EXTNAME"

    if bscale != 1.0:
        bintable.header.set("BSCALE", bscale, after=after_keyword)
        after_keyword = "BSCALE"

    if bzero != 0.0:
        bintable.header.set("BZERO", bzero, after=after_keyword)

    try:
        bitpix_comment = image_header.comments["BITPIX"]
    except (AttributeError, KeyError):
        bitpix_comment = "data type of original image"

    try:
        naxis_comment = image_header.comments["NAXIS"]
    except (AttributeError, KeyError):
        naxis_comment = "dimension of original image"

    # Set the label for the first column in the table

    bintable.header.set(
        "TTYPE1", "COMPRESSED_DATA", "label for field 1", after="TFIELDS"
    )

    # Set the data format for the first column.  It is dependent
    # on the requested compression type.

    if compression_type == "PLIO_1":
        tform1 = "1QI" if huge_hdu else "1PI"
    else:
        tform1 = "1QB" if huge_hdu else "1PB"

    bintable.header.set(
        "TFORM1",
        tform1,
        "data format of field: variable length array",
        after="TTYPE1",
    )

    # Create the first column for the table.  This column holds the
    # compressed data.
    col1 = Column(name=bintable.header["TTYPE1"], format=tform1)

    # Create the additional columns required for floating point
    # data and calculate the width of the output table.

    zbitpix = image_header["BITPIX"]

    if zbitpix < 0 and quantize_level != 0.0:
        # floating point image has 'COMPRESSED_DATA',
        # 'GZIP_COMPRESSED_DATA', 'ZSCALE', and 'ZZERO' columns (unless using
        # lossless compression, per CFITSIO)
        ncols = 4

        ttype2 = "GZIP_COMPRESSED_DATA"

        # The required format for the GZIP_COMPRESSED_DATA is actually
        # missing from the standard docs, but CFITSIO suggests it
        # should be 1PB, which is logical.
        tform2 = "1QB" if huge_hdu else "1PB"

        # Set up the second column for the table that will hold any
        # uncompressable data.
        bintable.header.set("TTYPE2", ttype2, "label for field 2", after="TFORM1")

        bintable.header.set(
            "TFORM2",
            tform2,
            "data format of field: variable length array",
            after="TTYPE2",
        )

        col2 = Column(name=ttype2, format=tform2)

        # Set up the third column for the table that will hold
        # the scale values for quantized data.
        bintable.header.set("TTYPE3", "ZSCALE", "label for field 3", after="TFORM2")
        bintable.header.set(
            "TFORM3", "1D", "data format of field: 8-byte DOUBLE", after="TTYPE3"
        )
        col3 = Column(name=bintable.header["TTYPE3"], format=bintable.header["TFORM3"])

        # Set up the fourth column for the table that will hold
        # the zero values for the quantized data.
        bintable.header.set("TTYPE4", "ZZERO", "label for field 4", after="TFORM3")
        bintable.header.set(
            "TFORM4", "1D", "data format of field: 8-byte DOUBLE", after="TTYPE4"
        )
        after = "TFORM4"
        col4 = Column(name=bintable.header["TTYPE4"], format=bintable.header["TFORM4"])

        # Create the ColDefs object for the table
        cols = ColDefs([col1, col2, col3, col4])
    else:
        # default table has just one 'COMPRESSED_DATA' column
        ncols = 1
        after = "TFORM1"

        # Create the ColDefs object for the table
        cols = ColDefs([col1])

    # Update the table header with the width of the table, the
    # number of fields in the table, the indicator for a compressed
    # image HDU, the data type of the image data and the number of
    # dimensions in the image data array.
    bintable.header.set("NAXIS1", cols.dtype.itemsize, "width of table in bytes")
    bintable.header.set(
        "TFIELDS", ncols, "number of fields in each row", after="GCOUNT"
    )
    bintable.header.set(
        "ZIMAGE", True, "extension contains compressed image", after=after
    )
    bintable.header.set("ZBITPIX", zbitpix, bitpix_comment, after="ZIMAGE")
    bintable.header.set("ZNAXIS", image_header["NAXIS"], naxis_comment, after="ZBITPIX")

    # Verify that any input tile size parameter is the appropriate
    # size to match the HDU's data.

    tile_shape = _validate_tile_shape(
        tile_shape=tile_shape,
        compression_type=compression_type,
        image_header=image_header,
    )

    # Set up locations for writing the next cards in the header.
    last_znaxis = "ZNAXIS"

    if image_header["NAXIS"] > 0:
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

        if naxis in image_header:
            bintable.header.set(
                znaxis, axis, image_header.comments[naxis], after=last_znaxis
            )
        else:
            bintable.header.set(
                znaxis, axis, "length of original image axis", after=last_znaxis
            )

        bintable.header.set(ztile, ts, "size of tiles to be compressed", after=after1)
        last_znaxis = znaxis
        after1 = ztile

    # Set the NAXIS2 header card in the table hdu to the number of
    # rows in the table.
    bintable.header.set("NAXIS2", nrows, "number of rows in table")

    # Set the compression parameters in the table header based on the compression type

    after_keyword = "ZCMPTYPE"
    idx = 1

    if compression_type == "RICE_1":
        bintable.header.set(
            "ZNAME1", "BLOCKSIZE", "compression block size", after=after_keyword
        )
        bintable.header.set(
            "ZVAL1", DEFAULT_BLOCK_SIZE, "pixels per block", after="ZNAME1"
        )

        bintable.header.set(
            "ZNAME2", "BYTEPIX", "bytes per pixel (1, 2, 4, or 8)", after="ZVAL1"
        )

        if bintable.header["ZBITPIX"] == 8:
            bytepix = 1
        elif bintable.header["ZBITPIX"] == 16:
            bytepix = 2
        else:
            bytepix = DEFAULT_BYTE_PIX

        bintable.header.set(
            "ZVAL2", bytepix, "bytes per pixel (1, 2, 4, or 8)", after="ZNAME2"
        )
        after_keyword = "ZVAL2"
        idx = 3
    elif compression_type == "HCOMPRESS_1":
        bintable.header.set(
            "ZNAME1", "SCALE", "HCOMPRESS scale factor", after=after_keyword
        )
        bintable.header.set(
            "ZVAL1", hcomp_scale, "HCOMPRESS scale factor", after="ZNAME1"
        )
        bintable.header.set(
            "ZNAME2", "SMOOTH", "HCOMPRESS smooth option", after="ZVAL1"
        )
        bintable.header.set(
            "ZVAL2", hcomp_smooth, "HCOMPRESS smooth option", after="ZNAME2"
        )
        after_keyword = "ZVAL2"
        idx = 3

    if image_header["BITPIX"] < 0:  # floating point image
        bintable.header.set(
            "ZNAME" + str(idx),
            "NOISEBIT",
            "floating point quantization level",
            after=after_keyword,
        )
        bintable.header.set(
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

            bintable.header.set(
                "ZQUANTIZ",
                QUANTIZE_METHOD_NAMES[quantize_method],
                zquantiz_comment,
                after="ZVAL" + str(idx),
            )
        else:
            # If the ZQUANTIZ keyword is missing the default is to assume
            # no dithering, rather than whatever DEFAULT_QUANTIZE_METHOD
            # is set to
            quantize_method = bintable.header.get("ZQUANTIZ", NO_DITHER)

            if isinstance(quantize_method, str):
                for k, v in QUANTIZE_METHOD_NAMES.items():
                    if v.upper() == quantize_method:
                        quantize_method = k
                        break
                else:
                    quantize_method = NO_DITHER

        if quantize_method == NO_DITHER:
            if "ZDITHER0" in bintable.header:
                # If dithering isn't being used then there's no reason to
                # keep the ZDITHER0 keyword
                del bintable.header["ZDITHER0"]
        else:
            if dither_seed:
                dither_seed = generate_dither_seed(dither_seed)
            elif "ZDITHER0" in bintable.header:
                dither_seed = bintable.header["ZDITHER0"]
            else:
                dither_seed = generate_dither_seed(DEFAULT_DITHER_SEED)

            bintable.header.set(
                "ZDITHER0",
                dither_seed,
                "dithering offset when quantizing floats",
                after="ZQUANTIZ",
            )

    if image_header:
        # Move SIMPLE card from the image header to the
        # table header as ZSIMPLE card.

        if "SIMPLE" in image_header:
            bintable.header.set(
                "ZSIMPLE",
                image_header["SIMPLE"],
                image_header.comments["SIMPLE"],
                before="ZBITPIX",
            )

        # Move EXTEND card from the image header to the
        # table header as ZEXTEND card.

        if "EXTEND" in image_header:
            bintable.header.set(
                "ZEXTEND", image_header["EXTEND"], image_header.comments["EXTEND"]
            )

        # Move BLOCKED card from the image header to the
        # table header as ZBLOCKED card.

        if "BLOCKED" in image_header:
            bintable.header.set(
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
            bintable.header.set(
                "ZTENSION",
                "IMAGE",
                image_header.comments["XTENSION"],
                before="ZBITPIX",
            )

        # Move PCOUNT and GCOUNT cards from image header to the table
        # header as ZPCOUNT and ZGCOUNT cards.

        if "PCOUNT" in image_header:
            bintable.header.set(
                "ZPCOUNT",
                image_header["PCOUNT"],
                image_header.comments["PCOUNT"],
                after=last_znaxis,
            )

        if "GCOUNT" in image_header:
            bintable.header.set(
                "ZGCOUNT",
                image_header["GCOUNT"],
                image_header.comments["GCOUNT"],
                after="ZPCOUNT",
            )

        # Move CHECKSUM and DATASUM cards from the image header to the
        # table header as XHECKSUM and XDATASUM cards.

        if "CHECKSUM" in image_header:
            bintable.header.set(
                "ZHECKSUM",
                image_header["CHECKSUM"],
                image_header.comments["CHECKSUM"],
            )

        if "DATASUM" in image_header:
            bintable.header.set(
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
        if "XTENSION" in image_header:
            bintable.header.set(
                "ZTENSION",
                "IMAGE",
                image_header.comments["XTENSION"],
                before="ZBITPIX",
            )

        # Move PCOUNT and GCOUNT cards from image header to the table
        # header as ZPCOUNT and ZGCOUNT cards.

        if "PCOUNT" in image_header:
            bintable.header.set(
                "ZPCOUNT",
                image_header["PCOUNT"],
                image_header.comments["PCOUNT"],
                after=last_znaxis,
            )

        if "GCOUNT" in image_header:
            bintable.header.set(
                "ZGCOUNT",
                image_header["GCOUNT"],
                image_header.comments["GCOUNT"],
                after="ZPCOUNT",
            )

    # When we have an image checksum we need to ensure that the same
    # number of blank cards exist in the table header as there were in
    # the image header.  This allows those blank cards to be carried
    # over to the image header when the hdu is uncompressed.

    if "ZHECKSUM" in bintable.header:
        required_blanks = image_header._countblanks()
        image_blanks = image_header._countblanks()
        table_blanks = bintable.header._countblanks()

        for _ in range(required_blanks - image_blanks):
            image_header.append()
            table_blanks += 1

        for _ in range(required_blanks - table_blanks):
            bintable.header.append()

    bintable.columns = cols

    # Add any keywords that are in the original header that are not already
    # FIXME: don't use keyword_remaps, instead define an actual list to check
    # including regular expressions for NAXIS and other similar keywords
    for card in image_header.cards:
        if card.keyword == "":
            continue  # TODO: investigate why some cards have empty keyword
        elif card.keyword == "COMMENT":
            bintable.header.add_comment(card.value)
        elif card.keyword == "HISTORY":
            bintable.header.add_history(card.value)
        elif (
            card.keyword not in CompImageHeader._keyword_remaps
            and card.keyword not in bintable.header
            and not card.keyword.startswith("NAXIS")
        ):
            bintable.header.append(card)

    # TODO: avoid writing the same comment multiple times

    return bintable
