.. _ccddata:


CCDData Class
=============

Getting Started
---------------

Getting Data In
+++++++++++++++

Creating a `~astropy.nddata.CCDData` object from any array-like data using
`astropy.nddata` is convenient:

    >>> import numpy as np
    >>> from astropy.nddata import CCDData
    >>> ccd = CCDData(np.arange(10), unit="adu")

Note that behind the scenes, this creates references to (not copies of) your
data when possible, so modifying the data in ``ccd`` will modify the
underlying data.

You are **required** to provide a unit for your data. The most frequently used
units for these objects are likely to be ``adu``, ``photon``, and ``electron``,
which can be set either by providing the string name of the unit (as in the
example above) or from unit objects:

    >>> from astropy import units as u
    >>> ccd_photon = CCDData([1, 2, 3], unit=u.photon)
    >>> ccd_electron = CCDData([1, 2, 3], unit="electron")

If you prefer *not* to use the unit functionality, then use the special unit
``u.dimensionless_unscaled`` when you create your `~astropy.nddata.CCDData`
images:

    >>> ccd_unitless = CCDData(np.zeros((10, 10)),
    ...                        unit=u.dimensionless_unscaled)

A `~astropy.nddata.CCDData` object can also be initialized from a FITS file:

    >>> ccd = CCDData.read('my_file.fits', unit="adu")  # doctest: +SKIP

If there is a unit in the FITS file (in the ``BUNIT`` keyword), that will be
used, but explicitly providing a unit in ``read`` will override any unit in the
FITS file.

There is no restriction at all on what the unit can be — any unit in
`astropy.units` or another that you create yourself will work.

In addition, the user can specify the extension in a FITS file to use:

    >>> ccd = CCDData.read('my_file.fits', hdu=1, unit="adu")  # doctest: +SKIP

If ``hdu`` is not specified, it will assume the data is in the primary
extension. If there is no data in the primary extension, the first extension
with image data will be used.

Metadata
++++++++

When initializing from a FITS file, the ``header`` property is initialized using
the header of the FITS file. Metadata is optional, and can be provided by any
dictionary or dict-like object:

    >>> ccd_simple = CCDData(np.arange(10), unit="adu")
    >>> my_meta = {'observer': 'Edwin Hubble', 'exposure': 30.0}
    >>> ccd_simple.header = my_meta  # or use ccd_simple.meta = my_meta

Whether the metadata is case-sensitive or not depends on how it is
initialized. A FITS header, for example, is not case-sensitive, but a Python
dictionary is.

Getting Data Out
++++++++++++++++

A `~astropy.nddata.CCDData` object behaves like a ``numpy`` array (masked if the
`~astropy.nddata.CCDData` mask is set) in expressions, and the underlying
data (ignoring any mask) is accessed through the ``data`` attribute:

    >>> ccd_masked = CCDData([1, 2, 3], unit="adu", mask=[0, 0, 1])
    >>> 2 * np.ones(3) * ccd_masked   # one return value will be masked
    masked_array(data=[2.0, 4.0, --],
                 mask=[False, False,  True],
           fill_value=1e+20)
    >>> 2 * np.ones(3) * ccd_masked.data   # ignores the mask  # doctest: +FLOAT_CMP
    array([2., 4., 6.])

You can force conversion to a ``numpy`` array with:

    >>> np.asarray(ccd_masked)
    array([1, 2, 3])
    >>> np.ma.array(ccd_masked.data, mask=ccd_masked.mask)
    masked_array(data=[1, 2, --],
                 mask=[False, False,  True],
           fill_value=999999)

A method for converting a `~astropy.nddata.CCDData` object to a FITS HDU list
is also available. It converts the metadata to a FITS header:

    >>> hdulist = ccd_masked.to_hdu()

You can also write directly to a FITS file:

    >>> ccd_masked.write('my_image.fits')

Masks and Flags
+++++++++++++++

Although it is not required when a `~astropy.nddata.CCDData` image is created,
you can also specify a mask and/or flags.

A mask is a boolean array the same size as the data in which a value of
``True`` indicates that a particular pixel should be masked (*i.e.*, not be
included in arithmetic operations or aggregation).

Flags are one or more additional arrays (of any type) whose shape matches the
shape of the data. One particularly useful type of flag is a bit planes; for
more details about bit planes and the functions ``astropy`` provides for
converting them to binary masks, see :ref:`bitmask_details`. For more details
on setting flags, see `~astropy.nddata.NDData`.

WCS
+++

The ``wcs`` attribute of a `~astropy.nddata.CCDData` object can be set two ways.

+ If the `~astropy.nddata.CCDData` object is created from a FITS file that has
  WCS keywords in the header, the ``wcs`` attribute is set to a
  `~astropy.wcs.WCS` object using the information in the FITS header.

+ The WCS can also be provided when the `~astropy.nddata.CCDData` object is
  constructed with the ``wcs`` argument.

Either way, the ``wcs`` attribute is kept up to date if the
`~astropy.nddata.CCDData` image is trimmed.

Uncertainty
-----------

You can set the uncertainty directly, either by creating a
`~astropy.nddata.StdDevUncertainty` object first:

    >>> data = np.random.normal(size=(10, 10), loc=1.0, scale=0.1)
    >>> ccd = CCDData(data, unit="electron")
    >>> from astropy.nddata.nduncertainty import StdDevUncertainty
    >>> uncertainty = 0.1 * ccd.data  # can be any array whose shape matches the data
    >>> my_uncertainty = StdDevUncertainty(uncertainty)
    >>> ccd.uncertainty = my_uncertainty

Or by providing a `~numpy.ndarray` with the same shape as the data:

    >>> ccd.uncertainty = 0.1 * ccd.data  # doctest: +ELLIPSIS
    INFO: array provided for uncertainty; assuming it is a StdDevUncertainty. [...]

In this case, the uncertainty is assumed to be
`~astropy.nddata.StdDevUncertainty`.

Two other uncertainty classes are available for which error propagation is
also supported: `~astropy.nddata.VarianceUncertainty` and
`~astropy.nddata.InverseVariance`. Using one of these three uncertainties is
required to enable error propagation in `~astropy.nddata.CCDData`.

If you want access to the underlying uncertainty, use its ``.array`` attribute:

    >>> ccd.uncertainty.array  # doctest: +ELLIPSIS
    array(...)

Arithmetic with Images
----------------------

Methods are provided to perform arithmetic operations with a
`~astropy.nddata.CCDData` image and a number, an ``astropy``
`~astropy.units.Quantity` (a number with units), or another
`~astropy.nddata.CCDData` image.

Using these methods propagates errors correctly (if the errors are
uncorrelated), takes care of any necessary unit conversions, and applies masks
appropriately. Note that the metadata of the result is *not* set if the
operation is between two `~astropy.nddata.CCDData` objects.

    >>> result = ccd.multiply(0.2 * u.adu)
    >>> uncertainty_ratio = result.uncertainty.array[0, 0]/ccd.uncertainty.array[0, 0]
    >>> round(uncertainty_ratio, 5)   # doctest: +FLOAT_CMP
    0.2
    >>> result.unit
    Unit("adu electron")

.. note::
    The affiliated package `ccdproc <https://ccdproc.readthedocs.io>`_ provides
    functions for many common data reduction operations. Those functions try to
    construct a sensible header for the result and provide a mechanism for
    logging the action of the function in the header.


The arithmetic operators ``*``, ``/``, ``+``, and ``-`` are *not* overridden.

.. note::
   If two images have different WCS values, the ``wcs`` on the first
   `~astropy.nddata.CCDData` object will be used for the resultant object.
