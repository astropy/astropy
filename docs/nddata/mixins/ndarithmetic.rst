Arithmetic mixin
================

TODO: example(s)

Provided that the world coordinate system (WCS) and shape match, the
`~astropy.nddata.NDArithmetic` mixin provides methods
:meth:`~astropy.nddata.NDArithmetic.add`,
:meth:`~astropy.nddata.NDArithmetic.subtract`,
:meth:`~astropy.nddata.NDArithmetic.multiply` and
:meth:`~astropy.nddata.NDArithmetic.divide` methods is to allow the
combination of two data objects that have common WCS and shape and units
consistent with the operation performed, with consistent behavior for masks,
and with a framework to propagate uncertainties. These methods are intended
for use by sub-classes and functions that deal with more complex combinations.

Entries that are masked in either of the operands are also masked in the
result.

The :class:`~astropy.nddata.NDDataArithmetic` class is based on `~astropy.nddata.NDData` and includes arithmetic.

.. warning:: Uncertainty propagation is still experimental, and does not take
             into account correlated uncertainties.
