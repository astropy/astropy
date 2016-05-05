.. _nddata_io:

NDData File I/O
===============

Introduction
------------

`~astropy.nddata.NDDataRef` has two methods for reading and writing:

- :meth:`~astropy.nddata.NDIOMixin.read`
- :meth:`~astropy.nddata.NDIOMixin.write`

These use the astropy :ref:`io_registry`.

Currently only support for reading and writing ``.fits`` files is implemented.

Reading and writing FITS files
------------------------------

.. warning::
    By default everything that was written by
    :meth:`~astropy.nddata.NDIOMixin.write` should be readable using
    :meth:`~astropy.nddata.NDIOMixin.read`.

    additional requirements:

    - ``mask``, if set, needs to be a `numpy.ndarray` or convertible to one.
    - ``uncertainty.array``, if set, needs to be a `numpy.ndarray` or convertible to one.
    - ``wcs``, if set, needs to be a ``astropy.wcs.WCS`` object.

    The ``.fits`` reader and writer are new in astropy v.1.2 and may be changed
    in subsequent versions if the need arises.


To use these you have to include a ``format='simple_fits'`` parameter when
calling these methods::

    >>> import numpy as np
    >>> from astropy.nddata import NDDataRef

    >>> ndd = NDDataRef([1,2,3,4])
    >>> ndd.write('filename.fits', format='simple_fits')

and reading the file again

    >>> ndd2 = NDDataRef.read('filename.fits', format='simple_fits')
    >>> ndd2
    NDDataRef([1, 2, 3, 4])


Reading external FITS files
---------------------------

External files probably don't follow the conventions used by
:meth:`~astropy.nddata.NDIOMixin.write` so it may be needed to use additional
parameters while reading such a file.

Specify different extensions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default :meth:`~astropy.nddata.NDIOMixin.read` assumes:

- ``data`` is saved in the primary HDU: extension ``0``.
- ``meta`` is saved in the primary HDU: extension ``0``.
- ``mask`` is saved in an extension called ``"MASK"`` or not present.
- ``uncertainty`` is saved in an extension called ``"UNCERT"`` or not present.
- ``unit`` is saved as value in the ``meta`` with keyword ``"BUNIT"`` or not
  present.

If the specified extension for ``mask`` or ``uncertainty`` or the
keyword for the ``unit`` doesn't exist there will be **NO** warning or info
message.

To overwrite any of these defaults use:

- ``ext_data``, ``ext_meta``, ``ext_mask``, ``ext_uncert`` parameter if the
  attribute is saved in another extension. Specifying extensions in ``read``
  and ``write`` can be given as string or as number but **not** both.

- ``kw_unit`` if the keyword for the unit differs from the default.

For example the ``data`` is saved in an extension called ``"SCI"``::

    >>> ndd = NDDataRef.read('filename.fits', format='simple_fits', ext_data='sci') # doctest: +SKIP

.. warning::
    Reading a compressed not-array-like ``mask`` is not possible.

Special cases
^^^^^^^^^^^^^

Since uncertainties have to be wrapped in an appropriate class you might need
to do  an additional step afterwards::

    >>> ndd = NDDataRef.read('filename.fits', format='simple_fits') # doctest: +SKIP
    >>> # For example if it's a standard deviation uncertainty:
    >>> ndd.uncertainty = StdDevUncertainty(ndd.uncertainty.array) # doctest: +SKIP

Also there might be units that cannot be interpreted by `~astropy.units.Unit`.
In these cases you should set the parameter to ``None`` and manually add the
unit, if required, later::

    >>> ndd = NDDataRef.read('filename.fits', format='simple_fits', kw_unit=None) # doctest: +SKIP
    >>> ndd.meta # if you want to inspect the header # doctest: +SKIP
    >>> ndd.unit = 'adu' # doctest: +SKIP

In case you want to change the datatype of your data (maybe because the data
was saved in unsigned integer but you want floats) you can specify a ``dtype``
parameter::

    >>> ndd = NDDataRef.read('filename.fits', format='simple_fits', dtype=float) # doctest: +SKIP

this ``dtype`` will affect **only** the data. Other attributes like mask and
uncertainty will be unaffected. You can always manually alter their dtype using
the appropriate attribute setter::

    >>> ndd = NDDataRef([1,2,0], uncertainty=[1,2,3]) # doctest: +SKIP
    INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

    >>> # Change the data type of the uncertainty to float:
    >>> ndd.uncertainty = ndd.uncertainty.array.astype(float) # doctest: +SKIP
    INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

Additional parameters
^^^^^^^^^^^^^^^^^^^^^

The ``copy`` paramter (default ``False``) will be passed to
`~astropy.nddata.NDData` during initialization.

You can specify additional keywords that are passed to
:func:`astropy.io.fits.open`. Probably not all of these listed there might be
possible.

Writing FITS files
------------------

:meth:`~astropy.nddata.NDIOMixin.write` also supports some optional arguments
like ``ext_mask``, ``ext_uncert`` and ``kw_unit`` but generally it might not
be needed to use them if you don't need to process them using other software.

Additional parameters
^^^^^^^^^^^^^^^^^^^^^

Writing also supports giving parameters to
:meth:`astropy.io.fits.HDUList.writeto`. Especially ``clobber`` might be
helpful if replacing an existing file is desired::

    >>> ndd = NDDataRef([1,2,3,4]) # doctest: +SKIP
    >>> ndd.write('test.fits', format='simple_fits')  # doctest: +SKIP
    >>> ndd.data[1] = 100  # doctest: +SKIP
    >>> # Suppose you want to overwrite this file again use clobber=True
    >>> ndd.write('test.fits', format='simple_fits', clobber=True)  # doctest: +SKIP

Why simple?
-----------

FITS files come in a plethora of formats and with varying conventions. The
parameters for :meth:`~astropy.nddata.NDIOMixin.read` allow some flexibility
but these don't cover all cases. It may be easier to define a customized
reader and writer (for inspiration take a look at the source code of the
``"simple_fits"`` code in ``"astropy.io.fits.connect.py"``) to handle
incompatible formats.
