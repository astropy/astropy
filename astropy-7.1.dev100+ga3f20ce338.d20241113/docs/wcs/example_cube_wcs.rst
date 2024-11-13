.. _example_cube_wcs:

Second Example
^^^^^^^^^^^^^^

Another way of creating a WCS object is via the use of a Python
dictionary. This affords us more control over the ``NAXISn``
FITS header keyword which is otherwise automatically default to zero
as in the case of the First Example shown above.

.. literalinclude:: examples/cube_wcs.py
   :language: python
