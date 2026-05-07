"""Helpers for handling FITS logical (``'L'``) variable-length array data."""

import numpy as np


def _logical_to_fits_bytes(row):
    """Convert a logical-VLA row to its FITS L wire-format bytes.

    Returns an int8 array containing ord('T') / ord('F'). Bool inputs map
    to T/F; numeric inputs treat zero as False and non-zero as True.
    Bytes input (|S1, e.g. from reading with ``logical_as_bytes=True``)
    is viewed as int8 verbatim so that NULL (b'\\x00') survives.
    """
    arr = np.asarray(row)
    if arr.dtype.kind == "S":
        return arr.view(np.int8)
    if arr.dtype == bool:
        return np.where(arr, ord("T"), ord("F")).astype(np.int8)
    return np.where(arr == 0, ord("F"), ord("T")).astype(np.int8)


def _logical_vla_heap_has_null(raw_data, field, heap_offset):
    """Return True if the heap unambiguously contains NULL bytes.

    Used to decide whether to warn that NULL values would be silently
    converted to False when reading a logical VLA column without
    ``logical_as_bytes=True``. To avoid a misleading warning on an
    all-zero heap (which under the FITS standard means all-NULL but
    under the astropy <= 7.2.0 legacy encoding means all-False), this
    requires both a 0x00 byte AND at least one ord('T') / ord('F')
    byte. The legacy-detection heuristic in
    ``_detect_legacy_logical_vla_heap`` separately catches files with
    a 0x01 byte; the only case this still treats as ambiguous is a
    heap containing nothing but 0x00, where the read value (``False``)
    is correct under either interpretation.
    """
    has_null = False
    has_tf = False
    for idx in range(len(field)):
        offset = int(field[idx, 1]) + heap_offset
        count = int(field[idx, 0])
        if not count:
            continue
        chunk = raw_data[offset : offset + count]
        if not has_null and (chunk == 0).any():
            has_null = True
        if not has_tf and ((chunk == ord("T")) | (chunk == ord("F"))).any():
            has_tf = True
        if has_null and has_tf:
            return True
    return False


def _detect_legacy_logical_vla_heap(raw_data, field, heap_offset):
    """Heuristically detect a logical VLA column written by astropy <= 7.2.0.

    astropy <= 7.2.0 stored bool VLA values as 0x00/0x01 bytes instead of the
    FITS L wire format ord('T') (0x54) / ord('F') (0x46). A column is
    classified as legacy when its heap region contains the byte 0x01 and no
    bytes other than 0x00 / 0x01 (those values cannot occur in a correctly
    written L column, where only 0x00, 0x46, 0x54 are valid). Heaps
    containing only 0x00 are ambiguous but indistinguishable in either
    interpretation.
    """
    bufs = []
    for idx in range(len(field)):
        offset = int(field[idx, 1]) + heap_offset
        count = int(field[idx, 0])
        if count:
            bufs.append(raw_data[offset : offset + count].view(np.uint8))
    if not bufs:
        return False
    unique = np.unique(np.concatenate(bufs))
    return bool(np.all(unique <= 1)) and 1 in unique
