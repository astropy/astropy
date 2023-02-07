import numbers

__all__ = ["simplify_basic_index"]


def simplify_basic_index(basic_index, *, shape):
    """
    Given a Numpy basic index, return a tuple of integers and slice objects
    with no default values (`None`) if possible.

    If one of the dimensions has a slice and the step is negative and the stop
    value of the slice was originally `None`, the new stop value of the slice
    may still be set to `None`.

    For more information on valid basic indices, see
    https://numpy.org/doc/stable/user/basics.indexing.html#basic-indexing

    Parameters
    ----------
    basic_index
        A valid Numpy basic index
    shape
        The shape of the array being indexed
    """
    ndim = len(shape)

    if not isinstance(basic_index, (tuple, list)):  # We just have a single int
        basic_index = (basic_index,)

    new_index = list(basic_index)

    if Ellipsis in new_index:
        if new_index.count(Ellipsis) > 1:
            raise IndexError("an index can only have a single ellipsis ('...')")

        # Replace the Ellipsis with the correct number of slice(None)s
        e_ind = new_index.index(Ellipsis)
        new_index.remove(Ellipsis)
        n_e = ndim - len(new_index)
        for i in range(n_e):
            ind = e_ind + i
            new_index.insert(ind, slice(None))

    if len(new_index) > ndim:
        raise ValueError(
            f"The dimensionality of the basic index {basic_index} can not be greater "
            f"than the dimensionality ({ndim}) of the data."
        )

    for i in range(ndim):
        if i < len(new_index):
            slc = new_index[i]
            if isinstance(slc, slice):
                if slc.start is None:
                    if slc.step is not None and slc.step < 0:
                        start = shape[i]
                    else:
                        start = 0
                elif slc.start < 0:
                    start = shape[i] + slc.start
                else:
                    start = slc.start

                if slc.stop is None:
                    if slc.step is not None and slc.step < 0:
                        stop = None
                    else:
                        stop = shape[i]
                elif slc.stop < 0:
                    stop = shape[i] + slc.stop
                else:
                    stop = slc.stop

                if slc.step is None:
                    step = 1
                else:
                    step = slc.step

                start = min(start, shape[i] - 1)
                if stop is not None:
                    stop = min(stop, shape[i])

                new_index[i] = slice(start, stop, step)

            elif isinstance(slc, numbers.Integral):
                if slc < 0:
                    slc = shape[i] + slc

                new_index[i] = int(slc)

            else:
                raise RuntimeError(f"Unexpected index element: {slc}")

        else:
            new_index.append(slice(0, shape[i], 1))

    return tuple(new_index)
