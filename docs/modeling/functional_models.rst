.. _functional_models:

Functional Models
*****************

These are models that are mathematically motivated, generally as solutions to
mathematical problems.   See :doc:`physical_models` for physically motivated
models.

Simple Operations
-----------------

- :class:`~astropy.modeling.physical_models.Multiply` model multiples by a
  factor and propagates units if the factor is a :class:`~astropy.units.Quantity`.

- :class:`~astropy.modeling.physical_models.Scale` model multiples by a
  factor without changing the units of the result.

1D Models
---------

- :class:`~astropy.modeling.physical_models.Moffat1D`

2D Models
---------

- :class:`~astropy.modeling.physical_models.AiryDisk2D`
