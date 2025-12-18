Overview of Astropy File I/O
****************************

Astropy provides two main interfaces for reading and writing data:

- :ref:`High-level Unified I/O <io-unified>` interface that is designed to be consistent
  and easy to use. This allows working with different types of data such as
  :ref:`tables <astropy-table>`, :ref:`images <astropy_nddata>`, and even
  :ref:`cosmologies <cosmology_io>`.
- Low-level I/O sub-packages that are directly responsible for
  reading and writing data in specific formats such as :ref:`FITS <astropy-io-fits>`
  or :ref:`VOTable <astropy-io-votable>`.

In general we recommend starting with the high-level interface unless you have a
specific need for the low-level interface.

.. list-table:: Comparison of high-level and low-level interfaces
   :widths: 50 50
   :header-rows: 1

   * - High-level Unified I/O
     - Low-level I/O
   * - Use ``read()`` and ``write()`` class methods of
       output data class, e.g., ``data = QTable.read("data.fits")`` returns a
       `~astropy.table.QTable`.
     - Interfaces are specific to format, e.g., ``hdus = fits.open("data.fits")``
       returns an `~astropy.io.fits.HDUList`.
   * - Read and write entire file at once.
     - Support varies, e.g., :ref:`FITS <astropy-io-fits>` has memory-mapped
       read access.
   * - Automatically determine file format in common cases.
     - Specify format explicitly.
   * - Help documentation via class method, e.g., ``QTable.read.help("fits")``.
     - Help documentation varies, e.g., ``help(fits.open)`` or API docs.
