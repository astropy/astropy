.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

The `~astropy.nddata` package provides classes to represent images and other
gridded data, some essential functions for manipulating images, and the
infrastructure for package developers who wish to include support for the
image classes.

.. _astropy_nddata_getting_started:

Getting started
===============

NDData
------

The primary purpose of `~astropy.nddata.NDData` is to act as a *container* for
data, metadata, and other related information like a mask.

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional `numpy` array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd1 = NDData(array)

or something that can be converted to an `numpy.ndarray`::

    >>> ndd2 = NDData([1, 2, 3, 4])
    >>> ndd2
    NDData([1, 2, 3, 4])

and can be accessed again via the ``data`` attribute::

    >>> ndd2.data
    array([1, 2, 3, 4])

It also supports additional properties like a ``unit`` or ``mask`` for the
data, a ``wcs`` (world coordinate system) and ``uncertainty`` of the data and
additional ``meta`` attributes:

    >>> data = np.array([1,2,3,4])
    >>> mask = data > 2
    >>> unit = 'erg / s'
    >>> from astropy.nddata import StdDevUncertainty
    >>> uncertainty = StdDevUncertainty(np.sqrt(data)) # representing standard deviation
    >>> meta = {'object': 'fictional data.'}
    >>> from astropy.coordinates import SkyCoord
    >>> wcs = SkyCoord('00h42m44.3s', '+41d16m09s')
    >>> ndd = NDData(data, mask=mask, unit=unit, uncertainty=uncertainty,
    ...              meta=meta, wcs=wcs)
    >>> ndd
    NDData([1, 2, 3, 4])

The representation only displays the ``data``; the other attributes need to be
accessed directly, for example ``ndd.mask`` to access the mask.


NDDataRef
---------

Building upon this pure container `~astropy.nddata.NDDataRef` implements:

+ a ``read`` and ``write`` method to access astropy's unified file io interface.
+ simple arithmetics like addition, subtraction, division and multiplication.
+ slicing.

Instances are created in the same way::

    >>> from astropy.nddata import NDDataRef
    >>> ndd = NDDataRef(ndd)
    >>> ndd
    NDDataRef([1, 2, 3, 4])

But also support arithmetic (:ref:`nddata_arithmetic`) like addition::

    >>> import astropy.units as u
    >>> ndd2 = ndd.add([4, -3.5, 3, 2.5] * u.erg / u.s)
    >>> ndd2
    NDDataRef([ 5. , -1.5,  6. ,  6.5])

Because these operations have a wide range of options these are not available
using arithmetic operators like ``+``.

Slicing or indexing (:ref:`nddata_slicing`) is possible (issuing warnings if
some attribute cannot be sliced)::

    >>> ndd2[2:]  # discard the first two elements  # doctest: +FLOAT_CMP
    INFO: wcs cannot be sliced. [astropy.nddata.mixins.ndslicing]
    NDDataRef([6. , 6.5])
    >>> ndd2[1]   # get the second element  # doctest: +FLOAT_CMP
    INFO: wcs cannot be sliced. [astropy.nddata.mixins.ndslicing]
    NDDataRef(-1.5)


StdDevUncertainty
-----------------

Though the `~astropy.nddata` package supports any kind of gridded data, this
introduction will focus on the use of `~astropy.nddata` for two-dimensional
images. To get started, we'll construct a two-dimensional image with a few
sources, some Gaussian noise, and a "cosmic ray" which we will later mask out::

    >>> import numpy as np
    >>> from astropy.modeling.models import Gaussian2D
    >>> y, x = np.mgrid[0:500, 0:600]
    >>> data = (Gaussian2D(1, 150, 100, 20, 10, theta=0.5)(x, y) +
    ...         Gaussian2D(0.5, 400, 300, 8, 12, theta=1.2)(x,y) +
    ...         Gaussian2D(0.75, 250, 400, 5, 7, theta=0.23)(x,y) +
    ...         Gaussian2D(0.9, 525, 150, 3, 3)(x,y) +
    ...         Gaussian2D(0.6, 200, 225, 3, 3)(x,y))
    >>> data += 0.01 * np.random.randn(500, 600)
    >>> cosmic_ray_value = 0.997
    >>> data[100, 300:310] = cosmic_ray_value

This image has a large "galaxy" in the lower left and the "cosmic ray" is the
horizontal line in the lower middle of the image:

.. doctest-skip::

    >>> import matplotlib.pyplot as plt
    >>> plt.imshow(data, origin='lower')

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    y, x = np.mgrid[0:500, 0:600]
    data = (Gaussian2D(1, 150, 100, 20, 10, theta=0.5)(x, y) +
            Gaussian2D(0.5, 400, 300, 8, 12, theta=1.2)(x,y) +
            Gaussian2D(0.75, 250, 400, 5, 7, theta=0.23)(x,y) +
            Gaussian2D(0.9, 525, 150, 3, 3)(x,y) +
            Gaussian2D(0.6, 200, 225, 3, 3)(x,y))
    np.random.seed(123456)
    data += 0.01 * np.random.randn(500, 600)
    cosmic_ray_value = 0.997
    data[100, 300:310] = cosmic_ray_value
    plt.imshow(data, origin='lower')


The "cosmic ray" can be masked out, in this simple test image, like this::

    >>> mask = (data == cosmic_ray_value)

`~astropy.nddata.CCDData` class for images
------------------------------------------

The `~astropy.nddata.CCDData` object, like the other objects in this package,
can store the data, a mask, and metadata. The `~astropy.nddata.CCDData` object
requires that a unit be specified::

    >>> from astropy.nddata import CCDData
    >>> ccd = CCDData(data, mask=mask,
    ...               meta={'object': 'fake galaxy', 'filter': 'R'},
    ...               unit='adu')

Slicing
-------

Slicing the works the way you would expect, with the mask and, if present,
WCS, sliced appropriately also::

    >>> ccd2 = ccd[:200, :]
    >>> ccd2.data.shape
    (200, 600)
    >>> ccd2.mask.shape
    (200, 600)
    >>> # Show the mask in a region around the cosmic ray:
    >>> ccd2.mask[99:102, 299:311]
    array([[False, False, False, False, False, False, False, False, False,
            False, False, False],
           [False,  True,  True,  True,  True,  True,  True,  True,  True,
             True,  True, False],
           [False, False, False, False, False, False, False, False, False,
            False, False, False]]...)

For many applications it may be more convenient to use
`~astropy.nddata.Cutout2D`, described in `image_utilities`_.

Image arithmetic, including uncertainty
---------------------------------------

Methods are provided for basic arithmetic operations between images, including
propagation of uncertainties. Support for uncertainties is currently limited
to standard deviation. The example below creates an uncertainty that is simply
Poisson error; note that the masked version of the square root is used::

    >>> ccd.uncertainty = np.ma.sqrt(np.ma.abs(ccd.data))
    INFO: array provided for uncertainty; assuming it is a StdDevUncertainty. [astropy.nddata.ccddata]

If we make a copy of the image and add that to the original, the uncertainty
changes as expected::

    >>> ccd2 = ccd.copy()
    >>> added_ccds = ccd.add(ccd2, handle_meta='first_found')
    >>> added_ccds.uncertainty.array[0, 0] / ccd.uncertainty.array[0, 0] / np.sqrt(2) # doctest: +FLOAT_CMP
    0.99999999999999989

Reading and writing
-------------------

A `~astropy.nddata.CCDData` can be saved to a FITS file::

    >>> ccd.write('test_file.fits')

and can also be read in from a FITS file::

    >>> ccd2 = CCDData.read('test_file.fits')

Note the unit is stored in the ``BUNIT`` keyword in the header on saving, and is
read from the header if it is present.

.. _image_utilities:

Image utilities
---------------

Cutouts
^^^^^^^

Though slicing directly is one way to extract a subframe,
`~astropy.nddata.Cutout2D` provides more convenient access to cutouts from the
data. The example below pulls out the large "galaxy" in the lower left of the
image, with the center of the cutout at ``position``::

    >>> from astropy.nddata import Cutout2D
    >>> position = (149.7, 100.1)
    >>> size = (80, 100)     # pixels
    >>> cutout = Cutout2D(ccd, position, size)
    >>> plt.imshow(cutout.data, origin='lower') # doctest: +SKIP

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import CCDData
    from astropy.nddata import Cutout2D
    y, x = np.mgrid[0:500, 0:600]
    data = (Gaussian2D(1, 150, 100, 20, 10, theta=0.5)(x, y) +
            Gaussian2D(0.5, 400, 300, 8, 12, theta=1.2)(x,y) +
            Gaussian2D(0.75, 250, 400, 5, 7, theta=0.23)(x,y) +
            Gaussian2D(0.9, 525, 150, 3, 3)(x,y) +
            Gaussian2D(0.6, 200, 225, 3, 3)(x,y))
    np.random.seed(123456)
    data += 0.01 * np.random.randn(500, 600)
    cosmic_ray_value = 0.997
    data[100, 300:310] = cosmic_ray_value
    mask = (data == cosmic_ray_value)
    ccd = CCDData(data, mask=mask,
                  meta={'object': 'fake galaxy', 'filter': 'R'},
                  unit='adu')
    position = (149.7, 100.1)
    size = (80, 100)     # pixels
    cutout = Cutout2D(ccd, position, size)
    plt.imshow(cutout.data, origin='lower')

This cutout can also plot itself on the original image::

    >>> plt.imshow(ccd, origin='lower')  # doctest: +SKIP
    >>> cutout.plot_on_original(color='white') # doctest: +SKIP

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import CCDData, Cutout2D
    y, x = np.mgrid[0:500, 0:600]
    data = (Gaussian2D(1, 150, 100, 20, 10, theta=0.5)(x, y) +
            Gaussian2D(0.5, 400, 300, 8, 12, theta=1.2)(x,y) +
            Gaussian2D(0.75, 250, 400, 5, 7, theta=0.23)(x,y) +
            Gaussian2D(0.9, 525, 150, 3, 3)(x,y) +
            Gaussian2D(0.6, 200, 225, 3, 3)(x,y))
    np.random.seed(123456)
    data += 0.01 * np.random.randn(500, 600)
    cosmic_ray_value = 0.997
    data[100, 300:310] = cosmic_ray_value
    mask = (data == cosmic_ray_value)
    ccd = CCDData(data, mask=mask,
                  meta={'object': 'fake galaxy', 'filter': 'R'},
                  unit='adu')
    position = (149.7, 100.1)
    size = (80, 100)     # pixels
    cutout = Cutout2D(ccd, position, size)
    plt.imshow(ccd, origin='lower')
    cutout.plot_on_original(color='white')

The cutout also provides methods for find pixel coordinates in the original or
in the cutout; recall that ``position`` is the center of the cutout in the
original image::

    >>> position
    (149.7, 100.1)
    >>> cutout.to_cutout_position(position)  # doctest: +FLOAT_CMP
    (49.7, 40.099999999999994)
    >>> cutout.to_original_position((49.7, 40.099999999999994))  # doctest: +FLOAT_CMP
     (149.7, 100.1)

For more details, including constructing a cutout from world coordinates and
the options for handling cutouts that go beyond the bounds of the original
image, see :ref:`cutout_images`.

Image resizing
^^^^^^^^^^^^^^

The functions `~astropy.nddata.block_reduce` and
`~astropy.nddata.block_replicate` resize images. The example below reduces the
size of the image by a factor of 4. Note that the result is a `numpy.ndarray`;
the mask, metadata, etc are discarded:

.. doctest-requires:: skimage

    >>> from astropy.nddata import block_reduce, block_replicate
    >>> smaller = block_reduce(ccd, 4)
    >>> smaller
    array(...)
    >>> plt.imshow(smaller, origin='lower')  # doctest: +SKIP

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import block_reduce, block_replicate
    from astropy.nddata import CCDData, Cutout2D
    y, x = np.mgrid[0:500, 0:600]
    data = (Gaussian2D(1, 150, 100, 20, 10, theta=0.5)(x, y) +
            Gaussian2D(0.5, 400, 300, 8, 12, theta=1.2)(x,y) +
            Gaussian2D(0.75, 250, 400, 5, 7, theta=0.23)(x,y) +
            Gaussian2D(0.9, 525, 150, 3, 3)(x,y) +
            Gaussian2D(0.6, 200, 225, 3, 3)(x,y))
    np.random.seed(123456)
    data += 0.01 * np.random.randn(500, 600)
    cosmic_ray_value = 0.997
    data[100, 300:310] = cosmic_ray_value
    mask = (data == cosmic_ray_value)
    ccd = CCDData(data, mask=mask,
                  meta={'object': 'fake galaxy', 'filter': 'R'},
                  unit='adu')
    smaller = block_reduce(ccd.data, 4)
    plt.imshow(smaller, origin='lower')

By default, both `~astropy.nddata.block_reduce` and
`~astropy.nddata.block_replicate` conserve flux.

Other image classes
-------------------


There are two less restrictive classes, `~astropy.nddata.NDDataArray` and
`~astropy.nddata.NDDataRef`, that can be used to hold image data. They are
primarily of interest to those who may want to create their own image class by
subclassing from one of the classes in the `~astropy.nddata` package. The main
differences between them are:

+ `~astropy.nddata.NDDataRef` can be sliced and has methods for basic
  arithmetic operations, but the user needs to use one of the uncertainty
  classes to define an uncertainty. See :ref:`NDDataRef` for more detail.
  Most of its properties must be set when the object is created because they
  are not mutable.
+ `~astropy.nddata.NDDataArray` extends `~astropy.nddata.NDDataRef` by adding
  the methods necessary to all it to behave like a numpy array in expressions
  and adds setters for several properties. It lacks the ability to
  automatically recognize and read data from FITS files and does not attempt
  to automatically set the WCS property.
+ `~astropy.nddata.CCDData` extends `~astropy.nddata.NDDataArray` by setting
  up a default uncertainty class, sets up straightforward read/write to FITS
  files, automatically sets up a WCS property.


More general gridded data class
-------------------------------

There are two additional classes in the ``nddata`` package that are of
interest primarily to people that either need a custom image class that goes
beyond the classes discussed so far or who are working with gridded data that
is not an image.

+ `~astropy.nddata.NDData` is a container class for holding general gridded
  data. It includes a handful of basic attributes, but no slicing or arithmetic.
  More information about this class is in :ref:`nddata_details`.
+ `~astropy.nddata.NDDataBase` is an abstract base class that developers of new
  gridded data classes can subclass to declare that the new class follows the
  `~astropy.nddata.NDData` interface. More details are in
  :ref:`nddata_subclassing`.

Additional examples
===================

The list of packages below that use the ``nddata`` framework is intended to be
useful to either people writing their own image classes or for those looking
for an image class that goes beyond what `~astropy.nddata.CCDData` does.

+ The `SunPy project <http://sunpy.org/>`_ uses `~astropy.nddata.NDData` as the
  foundation for its
  `Map classes <http://docs.sunpy.org/en/stable/code_ref/map.html>`_.
+ The class `~astropy.nddata.NDDataRef` is used in
  `specutils <http://specutils.readthedocs.io/en/latest/>`_ as the basis for
  `Spectrum1D <http://specutils.readthedocs.io/en/latest/api/specutils.Spectrum1D.html>`_, which adds several methods useful for
  spectra.
+ The package `ndmapper <http://ndmapper.readthedocs.io/en/latest/>`_, which
  makes it easy to build reduction pipelines for optical data, uses
  `~astropy.nddata.NDDataArray` as its image object.
+ The package `ccdproc <http://ccdproc.readthedocs.io/en/latest/>`_ uses the
  `~astropy.nddata.CCDData` class throughout for implementing optical/IR image
  reduction.

Using ``nddata``
================

.. toctree::
   :maxdepth: 2

   ccddata.rst
   utils.rst
   decorator.rst
   nddata.rst
   mixins/index.rst
   subclassing.rst

.. note that if this section gets too long, it should be moved to a separate doc page - see the top of performance.inc for the instructions on how to do that
.. include:: performance.inc

Reference/API
=============

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:

.. automodapi:: astropy.nddata.utils
    :no-inheritance-diagram:

.. _APE 7: https://github.com/astropy/astropy-APEs/blob/master/APE7.rst
