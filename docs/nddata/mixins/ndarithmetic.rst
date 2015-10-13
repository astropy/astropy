Arithmetic mixin
================

Overview
--------

The `~astropy.nddata.NDArithmeticMixin` adds methods for performing basic
operations on `~astropy.nddata.NDData` objects:
:meth:`~astropy.nddata.NDArithmeticMixin.add`,
:meth:`~astropy.nddata.NDArithmeticMixin.subtract`,
:meth:`~astropy.nddata.NDArithmeticMixin.multiply` and
:meth:`~astropy.nddata.NDArithmeticMixin.divide`.

The operations are permitted only if the two operands have the same WCS and
shape and the units, if any, consistent with the operation performed.  The
result is masked at a particular grid point if either of the operands is
masked at that point.

The operations include a framework to propagate uncertainties that are based
on the classes `~astropy.nddata.NDUncertainty`.

.. warning:: Uncertainty propagation is still experimental, and does not take
             into account correlated uncertainties.

Usage
-----

As with other mixins, using the arithmetic mixin starts with defining your own
class::

    >>> from astropy.nddata import NDData, NDArithmeticMixin
    >>> class MyNDDataArithmetic(NDArithmeticMixin, NDData): pass

Then, create instances of this new object with your data the way you would
with a plain `~astropy.nddata.NDData` object::

    >>> ndd1 = MyNDDataArithmetic([1, 2])
    >>> ndd2 = MyNDDataArithmetic([3, 4])

The mixin adds methods on these instances for combining them::

    >>> ndd1.add(ndd2)
    MyNDDataArithmetic([4, 6])
    >>> ndd2.multiply(ndd1)
    MyNDDataArithmetic([3, 8])

One important note: the order you list the mixins and `~astropy.nddata.NDData`
matters; the base   class, `~astropy.nddata.NDData` should be on the far
right.
