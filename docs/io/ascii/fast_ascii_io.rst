.. include:: references.txt

.. _fast_ascii_io:

Fast ASCII I/O
--------------

While :mod:`astropy.io.ascii` was designed with flexibility and extensibility
in mind, there is also a less flexible but significantly faster Cython/C engine for
reading and writing ASCII files. By default, |read| and |write| will attempt to
use this engine when dealing with compatible formats. The following formats
are currently compatible with the fast engine:

 * ``basic``
 * ``commented_header``
 * ``csv``
 * ``no_header``
 * ``rdb``
 * ``tab``

The fast engine can also be enabled through the format parameter by prefixing
a compatible format with "fast" and then an underscore. In this case, |read|
will not fall back on an ordinary reader if fast reading fails.
For example::

   >>> from astropy.table import Table
   >>> t = ascii.read('file.csv', format='fast_csv')  # doctest: +SKIP
   >>> t.write('output.csv', format='ascii.fast_csv')  # doctest: +SKIP

To disable the fast engine, specify ``fast_reader=False`` or
``fast_writer=False``. For example::

   >>> t = ascii.read('file.csv', format='csv', fast_reader=False) # doctest: +SKIP
   >>> t.write('file.csv', format='csv', fast_writer=False) # doctest: +SKIP

.. Note:: Guessing and Fast reading

   By default |read| will try to guess the format of in the input data by successively
   trying different formats until one succeeds ([reference the guessing section]).
   For the default ``'ascii'`` format this means that a number of pure Python readers
   with no fast implementation will be tried before getting to the fast readers.

   **For optimum performance**, turn off guessing entirely (``guess=False``) or
   narrow down the format options as much as possible by specifying the format
   (e.g. ``format='csv'``) and/or other options such as the delimiter.

Reading
^^^^^^^
Since the fast engine is not part of the ordinary :mod:`astropy.io.ascii`
infrastructure, fast readers raise an error when passed certain
parameters which are not implemented in the fast reader
infrastructure. In this case |read| will fall back on the ordinary reader.
These parameters are:

 * Negative ``header_start`` (except for commented-header format)
 * Negative ``data_start``
 * ``data_start=None``
 * ``comment`` string not of length 1
 * ``delimiter`` string not of length 1
 * ``quotechar`` string not of length 1
 * ``converters``
 * ``Outputter``
 * ``Inputter``
 * ``data_Splitter``
 * ``header_Splitter``

Parallel and fast conversion options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In addition to ``True`` and ``False``, the parameter ``fast_reader`` can also
be a dict specifying one or both of two additional parameters, ``parallel`` and
``use_fast_converter``. For example::

   >>> ascii.read('data.txt', format='basic', fast_reader={'parallel': True, 'use_fast_converter': True}) # doctest: +SKIP

These options allow for even faster table reading when enabled, but both are
disabled by default because they come with some caveats.

The ``parallel`` parameter can be used to enable multiprocessing via
the ``multiprocessing`` module, and can either be set to a number (the number
of processes to use) or ``True``, in which case the number of processes will be
``multiprocessing.cpu_count()``.   Note that this can cause issues within the
IPython Notebook and so enabling multiprocessing in this context is discouraged.

Setting ``use_fast_converter`` to be ``True`` enables a faster but
slightly imprecise conversion method for floating-point values, as described below.

Writing
^^^^^^^
The fast engine supports the same functionality as the ordinary writing engine
and is generally about 2 to 4 times faster than the ordinary engine. An IPython
notebook testing the relative performance of the fast writer against the
ordinary writing system and the data analysis library `Pandas
<http://pandas.pydata.org/>`__ is available `here <http://nbviewer.ipython.org/github/astropy/astropy-notebooks/blob/master/io/ascii/ascii_write_bench.ipynb>`__.
The speed advantage of the faster engine is greatest for integer data and least
for floating-point data; the fast engine is around 3.6 times faster for a
sample file including a mixture of floating-point, integer, and text data.
Also note that stripping string values slows down the writing process, so
specifying ``strip_whitespace=False`` can improve performance.

Fast converter
^^^^^^^^^^^^^^
Input floating-point values should ideally be converted to the
nearest possible floating-point approximation; that is, the conversion
should be correct within half of the distance between the two closest
representable values, or 0.5 `ULP
<http://en.wikipedia.org/wiki/Unit_in_the_last_place>`__. The ordinary readers,
as well as the default fast reader, are guaranteed to convert floating-point
values within 0.5 ULP, but there is also a faster and less accurate
conversion method accessible via ``use_fast_converter``. If the input
data has less than about 15 significant figures, or if accuracy is relatively
unimportant, this converter might be the best option in
performance-critical scenarios.

`Here
<http://nbviewer.ipython.org/github/astropy/astropy-notebooks/blob/master/io/ascii/conversion_profile.ipynb>`__
is an IPython notebook analyzing the error of the fast converter, both in
decimal values and in ULP. For values with a reasonably small number of
significant figures, the fast converter is guaranteed to produce an optimal
conversion (within 0.5 ULP). Once the number of significant figures exceeds
the precision of 64-bit floating-point values, the fast converter is no
longer guaranteed to be within 0.5 ULP, but about 60% of values end up
within 0.5 ULP and about 90% within 1.0 ULP. Another notebook analyzing
the fast converter's behavior with extreme values (such as subnormals
and values out of the range of floats) is available `here
<http://nbviewer.ipython.org/github/astropy/astropy-notebooks/blob/master/io/ascii/test_converter.ipynb>`__.

Speed gains
^^^^^^^^^^^
The fast ASCII engine was designed based on the general parsing strategy
used in the `Pandas <http://pandas.pydata.org/>`__ data analysis library, so
its performance is generally comparable (although slightly slower by
default) to the Pandas ``read_csv`` method.
`Here
<http://nbviewer.ipython.org/github/astropy/astropy-notebooks/blob/master/io/ascii/ascii_read_bench.ipynb>`__
is an IPython notebook comparing the performance of the ordinary
:mod:`astropy.io.ascii` reader, the fast reader, the fast reader with the
fast converter enabled, numpy's ``genfromtxt``, and Pandas' ``read_csv``
for different kinds of table data in a basic space-delimited file.

In summary, ``genfromtxt`` and the ordinary :mod:`astropy.io.ascii` reader
are very similar in terms of speed, while ``read_csv`` is slightly faster
than the fast engine for integer and floating-point data; for pure
floating-point data, enabling the fast converter yields a speedup of about
50%. Also note that Pandas uses the exact same method as the fast
converter in AstroPy when converting floating-point data.

The difference in performance between the fast engine and Pandas for
text data depends on the extent to which data values are repeated, as
Pandas is almost twice as fast as the fast engine when every value is
identical and the reverse is true when values are randomized. This is
because the fast engine uses fixed-size numpy string arrays for
text data, while Pandas uses variable-size object arrays and uses an
underlying set to avoid copying repeated values.

Overall, the fast engine tends to be around 4 or 5 times faster than
the ordinary ASCII engine. If the input data is very large (generally
about 100,000 rows or greater), and particularly if the data doesn't
contain primarily integer data or repeated string values, specifying
``parallel`` as ``True`` can yield further performance gains. Although
IPython doesn't work well with ``multiprocessing``, there is a
`script <https://github.com/amras1/ascii-profiling/blob/master/parallel.py>`__
available for testing the performance of the fast engine in parallel,
and a sample result may be viewed `here
<http://amras1.github.io/ascii-profiling/>`__. This profile uses the
fast converter for both the serial and parallel AstroPy
readers.

Another point worth noting is that the fast engine uses memory mapping
if a filename is supplied as input. If you want to avoid this for whatever
reason, supply an open file object instead. However, this will generally
be less efficient from both a time and a memory perspective, as the entire
file input will have to be read at once.

