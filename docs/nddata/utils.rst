.. _nddata_utils:

astropy.nddata.utils
====================

Overview
--------

The `astropy.nddata.utils` module includes general utility functions
for array operations.

.. _cutout_images:

2D Cutout Images
----------------

Getting Started
^^^^^^^^^^^^^^^

The `~astropy.nddata.utils.Cutout2D` class can be used to create a
postage stamp cutout image from a 2D array.  If an optional
`~astropy.wcs.WCS` object is input to
`~astropy.nddata.utils.Cutout2D`, then the
`~astropy.nddata.utils.Cutout2D` object will contain an updated
`~astropy.wcs.WCS` corresponding to the cutout array.

First, let's simulate a single source on a 2D data array. If you would like to
simulate many sources, see :ref:`bounding-boxes`.

    >>> import numpy as np
    >>> from astropy.modeling.models import Gaussian2D
    >>> y, x = np.mgrid[0:500, 0:500]
    >>> data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)

Now, let's display the image:

.. doctest-skip::

    >>> import matplotlib.pyplot as plt
    >>> plt.imshow(data, origin='lower')

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    y, x = np.mgrid[0:500, 0:500]
    data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)
    plt.imshow(data, origin='lower')

Now let's create a cutout array for the single object in this image.
We create a cutout array centered at position ``(x, y) = (49.7,
100.1)`` with a shape of ``(ny, nx) = (40, 50)``::

    >>> from astropy.nddata import Cutout2D
    >>> from astropy import units as u
    >>> position = (49.7, 100.1)
    >>> shape = (40*u.pixel, 50*u.pixel)
    >>> cutout = Cutout2D(data, position, shape)

The cutout array is stored in the ``data`` attribute of the
`~astropy.nddata.utils.Cutout2D` instance.  If the ``copy`` keyword is
`False` (default), then ``cutout.data`` will be a view into the
original ``data`` array.  If ``copy=True``, then ``cutout.data`` will
hold a copy of the original ``data``.  Let's display the cutout
image:

.. doctest-skip::

    >>> plt.imshow(cutout.data, origin='lower')

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import Cutout2D
    y, x = np.mgrid[0:500, 0:500]
    data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)
    position = (49.7, 100.1)
    shape = (40, 50)
    cutout = Cutout2D(data, position, shape)
    plt.imshow(cutout.data, origin='lower')

The cutout object can plot its bounding box on the original data using
the :meth:`~astropy.nddata.utils.Cutout2D.plot_on_original` method:

.. doctest-skip::

    >>> plt.imshow(data, origin='lower')
    >>> cutout.plot_on_original(color='white')

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import Cutout2D
    y, x = np.mgrid[0:500, 0:500]
    data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)
    position = (49.7, 100.1)
    shape = (40, 50)
    cutout = Cutout2D(data, position, shape)
    plt.imshow(data, origin='lower')
    cutout.plot_on_original(color='white')

Many properties of the cutout array are also stored as attributes,
including::

    >>> # shape of the cutout array
    >>> print(cutout.shape)
    (40, 50)

    >>> # rounded pixel index of the input position
    >>> print(cutout.position_original)
    (50, 100)

    >>> # corresponding position in the cutout array
    >>> print(cutout.position_cutout)
    (25, 20)

    >>> # (non-rounded) input position in both the original and cutout arrays
    >>> print(cutout.input_position_original, cutout.input_position_cutout)    # doctest: +FLOAT_CMP
    ((49.7, 100.1), (24.700000000000003, 20.099999999999994))

    >>> # the origin pixel in both arrays
    >>> print(cutout.origin_original, cutout.origin_cutout)
    ((25, 80), (0, 0))

    >>> # tuple of slice objects for the original array
    >>> print(cutout.slices_original)
    (slice(80, 120, None), slice(25, 75, None))

    >>> # tuple of slice objects for the cutout array
    >>> print(cutout.slices_cutout)
    (slice(0, 40, None), slice(0, 50, None))

Cutouts don't have to be specified by their shape if they are square.
Let's create another cutout array centered at position ``(x, y) = (49.7,
100.1)``, but this time with a square cutout that is 55 pixels to a side::

    >>> side_length = 55*u.pixel
    >>> cutout2 = Cutout2D(data, position, side_length=side_length)

.. doctest-skip::

    >>> plt.imshow(cutout2.data, origin='lower')

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import Cutout2D
    y, x = np.mgrid[0:500, 0:500]
    data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)
    position = (49.7, 100.1)
    cutout = Cutout2D(data, position, side_length=55)
    plt.imshow(cutout.data, origin='lower')

There are also two `~astropy.nddata.utils.Cutout2D` methods to convert
pixel positions between the original and cutout arrays::

    >>> print(cutout.to_original_position((2, 1)))
    (27, 81)

    >>> print(cutout.to_cutout_position((27, 81)))
    (2, 1)


2D Cutout modes
^^^^^^^^^^^^^^^

There are three modes for creating cutout arrays, ``'trim'``,
``'partial'``, and ``'strict'``.  For the ``'partial'`` and ``'trim'``
modes, a partial overlap of the cutout array and the input ``data``
array is sufficient.  For the ``'strict'`` mode, the cutout array has
to be fully contained within the ``data`` array, otherwise an
`~astropy.nddata.utils.PartialOverlapError` is raised.   In all modes,
non-overlapping arrays will raise a
`~astropy.nddata.utils.NoOverlapError`.  In ``'partial'`` mode,
positions in the cutout array that do not overlap with the ``data``
array will be filled with ``fill_value``.  In ``'trim'`` mode only the
overlapping elements are returned, thus the resulting cutout array may
be smaller than the requested ``shape``.

The default uses ``mode='trim'``, which can result in cutout arrays
that are smaller than the requested ``shape``::

    >>> data2 = np.arange(20.).reshape(5, 4)
    >>> c1 = Cutout2D(data2, (0, 0), (3, 3), mode='trim')
    >>> print(c1.data)
    [[ 0.  1.]
     [ 4.  5.]]
    >>> print(c1.position_original, c1.position_cutout)
    ((0, 0), (0, 0))

With ``mode='partial'``, the cutout will never be trimmed.  Instead it
will be filled with ``fill_value`` (the default is ``numpy.nan``) if
the cutout is not fully contained in the data array::

    >>> c2 = Cutout2D(data2, (0, 0), (3, 3), mode='partial')
    >>> print(c2.data)
    [[ nan  nan  nan]
     [ nan   0.   1.]
     [ nan   4.   5.]]

Note that for the ``'partial'`` mode, the positions (and several other
attributes) are calculated for on the *valid* (non-filled) cutout
values::

    >>> print(c2.position_original, c2.position_cutout)
    ((0, 0), (1, 1))
    >>> print(c2.origin_original, c2.origin_cutout)
    ((0, 0), (1, 1))
    >>> print(c2.slices_original)
    (slice(0, 2, None), slice(0, 2, None))
    >>> print(c2.slices_cutout)
    (slice(1, 3, None), slice(1, 3, None))

Using ``mode='strict'`` will raise an exception if the cutout is not
fully contained in the data array:

.. doctest-skip::

    >>> c3 = Cutout2D(data2, (0, 0), (3, 3), mode='strict')
    PartialOverlapError: Arrays overlap only partially.


2D Cutout from a `~astropy.coordinates.SkyCoord` position
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The input ``position`` can also be specified as a
`~astropy.coordinates.SkyCoord`, in which case a `~astropy.wcs.WCS`
object must be input via the ``wcs`` keyword.

First, let's define a `~astropy.coordinates.SkyCoord` position and a
`~astropy.wcs.WCS` object for our data (usually this would come from
your FITS header)::

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy.wcs import WCS
    >>> position = SkyCoord('13h11m29.96s -01d19m18.7s', frame='icrs')
    >>> wcs = WCS(naxis=2)
    >>> rho = np.pi / 3.
    >>> scale = 0.05 / 3600.
    >>> wcs.wcs.cd = [[scale*np.cos(rho), -scale*np.sin(rho)],
    ...               [scale*np.sin(rho), scale*np.cos(rho)]]
    >>> wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    >>> wcs.wcs.crval = [position.ra.value, position.dec.value]
    >>> wcs.wcs.crpix = [50, 100]

Now let's create the cutout array using the
`~astropy.coordinates.SkyCoord` position and ``wcs`` object::

    >>> cutout = Cutout2D(data, position, (30, 40), wcs=wcs)
    >>> plt.imshow(cutout.data, origin='lower')   # doctest: +SKIP

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import Cutout2D
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    y, x = np.mgrid[0:500, 0:500]
    data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)
    position = SkyCoord('13h11m29.96s -01d19m18.7s', frame='icrs')
    shape = (40, 50)
    wcs = WCS(naxis=2)
    rho = np.pi / 3.
    scale = 0.05 / 3600.
    wcs.wcs.cd = [[scale*np.cos(rho), -scale*np.sin(rho)],
                  [scale*np.sin(rho), scale*np.cos(rho)]]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [position.ra.value, position.dec.value]
    wcs.wcs.crpix = [50, 100]
    cutout = Cutout2D(data, position, (30, 40), wcs=wcs)
    plt.imshow(cutout.data, origin='lower')

The ``wcs`` attribute of the `~astropy.nddata.utils.Cutout2D` object now
contains the propagated `~astropy.wcs.WCS` for the cutout array.
Let's find the sky coordinates for a given pixel in the cutout array.
Note that we need to use the ``cutout.wcs`` object for the cutout
positions::

    >>> from astropy.wcs.utils import pixel_to_skycoord
    >>> x_cutout, y_cutout = (5, 10)
    >>> pixel_to_skycoord(x_cutout, y_cutout, cutout.wcs)    # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (197.8747893, -1.32207626)>

We now find the corresponding pixel in the original ``data`` array and
its sky coordinates::

    >>> x_data, y_data = cutout.to_original_position((x_cutout, y_cutout))
    >>> pixel_to_skycoord(x_data, y_data, wcs)    # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (197.8747893, -1.32207626)>

As expected, the sky coordinates in the original ``data`` and the
cutout array agree.


Reference/API
=============

.. automodapi:: astropy.nddata.utils
    :no-inheritance-diagram:
