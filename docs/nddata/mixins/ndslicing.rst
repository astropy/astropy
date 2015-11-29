Slicing mixin
=============

The slicing mixin adds the ability to index an `~astropy.nddata.NDData` object
in a manner similar to a numpy array. Slicing returns a new
`~astropy.nddata.NDData` object with the shape indicated by the slice.

The first step in creating an `~astropy.nddata.NDData`-based object with
slicing is to create a new class::

    >>> from astropy.nddata import NDData, NDSlicingMixin
    >>> class MyNDDataWithSlicing(NDSlicingMixin, NDData): pass

Then, initialize the new object the same way you would initialize a plain
`~astropy.nddata.NDData` object::

    >>> sliceable_ndd = MyNDDataWithSlicing([1, 2, 3, 4])
    >>> sliceable_ndd[1:3]
    MyNDDataWithSlicing([2, 3])

One important note: the order you list the mixins and `~astropy.nddata.NDData`
matters; the base class, `~astropy.nddata.NDData` should be on the far right.
