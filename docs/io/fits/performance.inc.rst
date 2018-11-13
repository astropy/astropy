.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-io-fits-performance:

.. TODO: determine whether the following is quantitatively true, and either
.. uncomment or remove.

.. Performance Tips
.. ================
..
.. By default, :func:`astropy.io.fits.open` will open files using memory-mapping,
.. which means that the data is not necessarily read into memory until it is
.. needed. While memory-efficient, if memory is not a concern for you, you may
.. find that you can get better performance by turning memory mapping off, which
.. forces the data to be read into memory immediately:
..
.. .. doctest-skip::
..
..     >>> fits.open('example.fits', memmap=False)
