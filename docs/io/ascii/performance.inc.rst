.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-io-ascii-performance:

Performance Tips
================

By default, when trying to read a file the reader will guess the format, which
involves trying to read it with many different readers. For better performance
when dealing with large tables, it is recommended to specify the format and any
options explicitly, and turn off guessing as well.

Example
-------

..
  EXAMPLE START
  Performance Tips for Reading Large Tables with astropy.io.ascii

If you are reading a simple CSV file with a one-line header with column names,
the following::

    read('example.csv', format='basic', delimiter=',', guess=False)  # doctest: +SKIP

can be at least an order of magnitude faster than::

    read('example.csv')  # doctest: +SKIP

..
  EXAMPLE END
