"""
This module is for functions that do tricky things with Numpy arrays and dtypes
that are not normally supported in Numpy (but can work in limited cases
relevant to FITS) or that otherwise require workarounds.
"""


def realign_dtype(dtype, offsets):
    """
    Given a Numpy struct dtype object an a list of integer offsets, with one
    offset per field in the dtype, returns a new dtype where each field has the
    given offset.

    All offsets must be non-negative integers, but otherwise have no
    restrictions, and may overlap, per the usual rules for creating struct
    dtypes.  The new dtype will have an itemsize equal to the offset of the
    right-most field plus the width of that field.

    One restriction of this function is that it must not be used with object
    arrays--incorrect offsets may lead to invalid pointers in the arrays.
    However, this function is really only meant for use by astropy.io.fits and
    object arrays are not supported for FITS data anyhow.

    This function is used primarily to get around a shortcoming in Numpy that
    it is currently impossible to create dtypes with arbitrary offsets, *and*
    that have zero-width fields.  Both of these features are needed for full
    FITS support.  However, this will be fixed in a future version of Numpy at
    which point use of this hack can be deprecated.  See
    https://github.com/numpy/numpy/pull/6430
    """

    # Previously this was implemented in C, but then I realized that the C
    # version is not needed--the workaround is to use dtype.__setstate__
    # Note: There is a comment in the Numpy source code (see
    # https://github.com/numpy/numpy/blob/v1.10.1/numpy/core/src/multiarray/descriptor.c#L2226)
    # that this may be changed at some point.  But hopefully by then the fixes
    # in #6430 will be implemented, making this hack unnecessary to begin with.

    cls, args, state = dtype.__reduce__()

    names, fields = state[3:5]
    fields = fields.copy()

    itemsize = 0  # We will re-determine the itemsize based on the type
                  # of the field with the largest (offset + itemsize)

    if fields is None or len(offsets) != len(names):
        raise ValueError(
            "Dtype must be a structured dtype, and length of offsets list "
            "must be the same as the number of fields.")

    for name, offset in zip(names, offsets):
        field = fields[name]
        itemsize = max(itemsize, offset + field[0].itemsize)

        if offset != field[1]:
            fields[name] = (field[0], offset)

    new_typespec = '|V{0}'.format(itemsize)

    new_state = state[:4] + (fields, itemsize) + state[6:]

    new_dtype = cls(new_typespec, *args[1:])
    new_dtype.__setstate__(new_state)

    return new_dtype
