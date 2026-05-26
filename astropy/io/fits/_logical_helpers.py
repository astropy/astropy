"""Helpers for handling FITS logical (``'L'``) variable-length array data."""

import warnings

import numpy as np

from astropy.utils.exceptions import AstropyDeprecationWarning

_VALID_LOGICAL_BYTES = (ord("T"), ord("F"), 0)


def _validate_logical_input(row):
    """Validate user input for a FITS logical column.

    Accepts only bool arrays/sequences and |S1 byte arrays whose values
    are limited to ``b'T'``, ``b'F'``, and ``b'\\x00'`` (the FITS L
    wire-format values for True, False, and NULL).

    |S1 input with any other byte value raises ``ValueError``.

    Any other input type (int, float, object, ...) emits an
    ``AstropyDeprecationWarning`` and is left to the existing coercion
    path; in a future version such input will be rejected with
    ``TypeError``.
    """
    if isinstance(row, (bool, np.bool_)):
        return

    if isinstance(row, np.ndarray):
        arr = row
    else:
        try:
            arr = np.asarray(row)
        except Exception:
            _warn_logical_input(type(row).__name__)
            return

    if arr.dtype.kind == "S":
        if arr.dtype.itemsize > 1:
            raise ValueError(
                "FITS logical column byte input must have itemsize 1 "
                f"(|S1); got {arr.dtype.str}."
            )
        if arr.size:
            view = arr.view(np.uint8)
            allowed = (
                (view == _VALID_LOGICAL_BYTES[0])
                | (view == _VALID_LOGICAL_BYTES[1])
                | (view == _VALID_LOGICAL_BYTES[2])
            )
            if not bool(allowed.all()):
                bad = np.unique(view[~allowed])
                bad_repr = ", ".join(repr(bytes([int(b)])) for b in bad)
                raise ValueError(
                    "FITS logical column |S1 input contains invalid "
                    f"byte{'s' if bad.size > 1 else ''} {bad_repr}; "
                    "only b'T', b'F', and b'\\x00' are allowed."
                )
        return

    if arr.dtype == bool:
        return

    # Integer arrays (signed or unsigned, any width) are accepted
    # without warning iff every value is 0 or 1; those map unambiguously
    # to False / True. Internal callers that legitimately need to pass
    # arbitrary ``int8`` storage bytes (e.g. ``ColDefs._init_from_array``
    # reading from disk) use ``Column(..., _skip_validation=True)``.
    if arr.dtype.kind in ("i", "u"):
        if arr.size == 0 or bool(((arr == 0) | (arr == 1)).all()):
            return

    _warn_logical_input(arr.dtype.str)


def _warn_logical_input(typename):
    warnings.warn(
        "Passing input of type "
        f"{typename!r} to a FITS logical ('L') column is deprecated "
        "and will raise in a future version. Pass a bool array, or "
        "a |S1 byte array containing only b'T', b'F', or b'\\x00'.",
        AstropyDeprecationWarning,
        stacklevel=3,
    )


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
