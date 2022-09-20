# cython: language_level=3
from collections import OrderedDict


cdef Py_ssize_t BLOCK_SIZE = 2880  # the FITS block size
cdef Py_ssize_t CARD_LENGTH = 80
cdef str VALUE_INDICATOR = '= '  # The standard FITS value indicator
cdef str END_CARD = 'END' + ' ' * 77


def parse_header(fileobj):
    """Fast (and incomplete) parser for FITS headers.

    This parser only reads the standard 8 character keywords, and ignores the
    CONTINUE, COMMENT, HISTORY and HIERARCH cards. The goal is to find quickly
    the structural keywords needed to build the HDU objects.

    The implementation is straightforward: first iterate on the 2880-bytes
    blocks, then iterate on the 80-bytes cards, find the value separator, and
    store the parsed (keyword, card image) in a dictionary.

    """

    cards = OrderedDict()
    cdef list read_blocks = []
    cdef int found_end = 0
    cdef bytes block
    cdef str header_str, block_str, card_image, keyword
    cdef Py_ssize_t idx, end_idx, sep_idx

    while found_end == 0:
        # iterate on blocks
        block = fileobj.read(BLOCK_SIZE)
        if not block or len(block) < BLOCK_SIZE:
            # header looks incorrect, raising exception to fall back to
            # the full Header parsing
            raise Exception

        block_str = block.decode('ascii')
        read_blocks.append(block_str)
        idx = 0
        while idx < BLOCK_SIZE:
            # iterate on cards
            end_idx = idx + CARD_LENGTH
            card_image = block_str[idx:end_idx]
            idx = end_idx

            # We are interested only in standard keyword, so we skip
            # other cards, e.g. CONTINUE, HIERARCH, COMMENT.
            if card_image[8:10] == VALUE_INDICATOR:
                # ok, found standard keyword
                keyword = card_image[:8].strip()
                cards[keyword.upper()] = card_image
            else:
                sep_idx = card_image.find(VALUE_INDICATOR, 0, 8)
                if sep_idx > 0:
                    keyword = card_image[:sep_idx]
                    cards[keyword.upper()] = card_image
                elif card_image == END_CARD:
                    found_end = 1
                    break

    # we keep the full header string as it may be needed later to
    # create a Header object
    header_str = ''.join(read_blocks)
    return header_str, cards
