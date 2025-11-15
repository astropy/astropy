Adding a structured numpy array to a Table as a single column now emits an
``AstropyDeprecationWarning``. In astropy 5.2 and later, structured arrays
will be added as ``Column`` objects by default. To use ``NdarrayMixin``
explicitly, wrap the array: ``data.view(NdarrayMixin)``. To use ``Column``,
wrap the array: ``Column(data)``.
