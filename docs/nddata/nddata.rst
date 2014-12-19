NDData overview
===============

Classes
-------

`~astropy.nddata` provides both data and uncertainty classes. The plain
`~astropy.nddata.NDData` classes is a container for data in array form and its
associated metadata. `~astropy.nddata.NDDataArithmetic` adds support for
arithmetic operations on `~astropy.nddata.NDData` objects that supports
propagation of uncertainties. Those uncertainties are represented by subclasses
of `~astropy.nddata.NDUncertainty`; those classes are responsible for the
actual propagation of uncertainties and will in most cases need to be defined
by the user.

Initializing
------------

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional Numpy array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd = NDData(array)

Note that the data in ``ndd`` is a reference to the original ``array``, so
changing the data in ``ndd`` will change the corresponding data in ``array``
in most circumstances.

An `~astropy.nddata.NDData` object can also be instantiated by passing it an
`~astropy.nddata.NDData` object:

    >>> ndd1 = NDData(array)
    >>> ndd2 = NDData(ndd1)

As above, the data in``ndd2`` is a reference to the data in ``ndd1``, so
changes to one will affect the other.

If an `~astropy.nddata.NDData` object is initialized with a masked numpy array
then the mask attribute, described in the next section, is automatically set
to the mask of the numpy array.

Mask
----

Values can be masked using the ``mask`` attribute, which should be a boolean
Numpy array with the same dimensions as the data, e.g.::

     >>> ndd_masked = NDData(ndd, mask = ndd.data > 0.9)
     INFO: Overwriting NDData's current mask with specified mask [astropy.nddata.nddata]

or by initializing with a masked array. A mask value of `True` indicates a
value that should be ignored, while a mask value of `False` indicates a valid
value.

Uncertainties
-------------

`~astropy.nddata.NDData` objects have an ``uncertainty`` attribute that can be
used to set the uncertainty on the data values. The ``uncertainty`` must have
an attribute ``uncertainty_type`` which is a string.

While not a requirement, the following ``uncertainty_type`` strings
are strongly recommended for common ways of specifying normal
distributions:

+ ``"std"``: if ``uncertainty`` stores the standard deviation/sigma
  (either a single value or on a per-pixel basis).
+ ``"var"``: if ``uncertainty`` stores the variance (either a single
  value or on a per-pixel basis).
+ ``"ivar"``: if ``uncertainty`` stores the inverse variance (either a
  single value or on a per-pixel basis).


.. note:: For information on creating your own uncertainty classes,
          see :doc:`subclassing`.

Meta-data
---------

The :class:`~astropy.nddata.NDData` class includes a ``meta`` attribute
that defaults to an empty dictionary, and can be used to set overall meta-data
for the dataset::

    ndd.meta['exposure_time'] = 340.
    ndd.meta['filter'] = 'J'

Elements of the meta-data dictionary can be set to any valid Python object::

    ndd.meta['history'] = ['calibrated', 'aligned', 'flat-fielded']

Converting to Numpy arrays
--------------------------

`~astropy.nddata.NDData` objects can also be easily converted to
numpy arrays::

    >>> import numpy as np
    >>> arr = np.array(ndd)
    >>> np.all(arr.data == mydataarray)  # doctest: +SKIP
    True

TODO: IS BELOW REALLY TRUE?

If a ``mask`` is defined, this will result in a `~numpy.ma.MaskedArray`, so
in all cases a useable `numpy.ndarray` or subclass will result. This allows
straightforward plotting of `~astropy.nddata.NDData` objects with 1-
and 2-dimensional datasets using Matplotlib::

    >>> from matplotlib import pyplot as plt  # doctest: +SKIP
    >>> plt.plot(ndd)  # doctest: +SKIP

This works because the Matplotlib plotting functions automatically convert
their inputs using `numpy.array`.
