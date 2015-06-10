.. _nddata_utils:

astropy.nddata.utils
====================

Overview
--------

The `astropy.nddata.utils` module includes general utility functions
for array operations.


Cutout Images
-------------

Getting Started
^^^^^^^^^^^^^^^

The `~astropy.nddata.utils.Cutout` class can be used to create a
postage stamp cutout image from a 2D array.  If an optional
`~astropy.wcs.WCS` object is input to `~astropy.nddata.utils.Cutout`,
then the `~astropy.nddata.utils.Cutout` object will contain an updated
`~astropy.wcs.WCS` corresponding to the cutout array.

First, let's create a simple 2D data array:

.. plot::
    :include-source:

    import numpy as np
    from astropy.modeling.models import Gaussian2D
    y, x = np.mgrid[0:500, 0:500]
    data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)
    plt.imshow(data, origin='lower')

Now let's create a cutout array for the single object in this image.
We create a cutout array centered at position ``(100.1, 49.7)`` (``(y,
x)``) with a shape of ``(40, 50)`` (``(ny, nx)``)::

    >>> from astropy.nddata import Cutout
    >>> position = (100.1, 49.7)
    >>> shape = (40, 50)
    >>> cutout = Cutout(data, position, shape)

The cutout array is stored in the ``data`` attribute of the
`~astropy.nddata.utils.Cutout` instance::

    >>> plt.imshow(cutout.data, origin='lower')

.. plot::

    import numpy as np
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import Cutout
    y, x = np.mgrid[0:500, 0:500]
    data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)
    position = (100.1, 49.7)
    shape = (40, 50)
    cutout = Cutout(data, position, shape)
    plt.imshow(cutout.data, origin='lower')

Many other properties of the cutout array are also stored as
attributes, including::

    >>> # shape of the cutout array
    >>> print(cutout.shape)
    (40, 50)

    >>> # rounded pixel index of the input position
    >>> print(cutout.position_large)
    (100, 50)

    >>> # corresponding position in the cutout array
    >>> print(cutout.position_small)
    (20, 50)

    >>> # (non-rounded) input position in both the large and cutout arrays
    >>> print(cutout.input_position_large, cutout.input_position_small)
    ((100.1, 49.7), (20.099999999999994, 24.700000000000003))

    >>> # the origin pixel in both arrays
    >>> print(cutout.origin_large, cutout.origin_small)
    ((80, 25), (0, 0))

    >>> # tuple of slice objects for the large array
    >>> print(cutout.slices_large)
    ((slice(80, 120, None), slice(25, 75, None))

    >>> # tuple of slice objects for the cutout array
    >>> print(cutout.slices_small)
    (slice(0, 40, None), slice(0, 50, None))

There are also two `~astropy.nddata.utils.Cutout` methods to convert
pixel positions between the large and cutout arrays::

    >>> print(cutout.to_large((1, 2)))
    (81, 27)

    >>> print(cutout.from_large((81, 27)))
    (1, 2)


Cutout modes
^^^^^^^^^^^^

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
    >>> c1 = Cutout(data2, (0, 0), (3, 3), mode='trim')
    >>> print(c1.data)
    [[ 0.  1.]
     [ 4.  5.]]
    >>> print(c1.position_large, c1.position_small)
    ((0, 0), (0, 0))

With ``mode='partial'``, the cutout will never be trimmed.  Instead it
will be filled with ``fill_value`` (the default is ``numpy.nan``) if
the cutout is not fully contained in the data array::

    >>> c2 = Cutout(data2, (0, 0), (3, 3), mode='partial')
    >>> print(c2.data)
    [[ nan  nan  nan]
     [ nan   0.   1.]
     [ nan   4.   5.]]

Note that for the ``'partial'`` mode, the positions (and several other
attributes) are calculated for on the *valid* (non-filled) cutout
values::

    >>> print(c2.position_large, c2.position_small)
    ((0, 0), (1, 1))
    >>> print(c2.origin_large, c2.origin_small)
    ((0, 0), (1, 1))
    >>> print(c2.slices_large)
    (slice(0, 2, None), slice(0, 2, None))
    >>> print(c2.slices_small)
    (slice(1, 3, None), slice(1, 3, None))

Using ``mode='strict'`` will raise an exception if the cutout is not
fully contained in the data array::

    >>> c3 = Cutout(data2, (0, 0), (3, 3), mode='strict')
    PartialOverlapError: Arrays overlap only partially.


Cutout from a `~astropy.coordinates.SkyCoord` position
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    >>> wcs.wcs.crpix = [100, 50]

Now let's create the cutout array using the
`~astropy.coordinates.SkyCoord` position and ``wcs`` object::

    >>> cutout = Cutout(data, position, (30, 40), wcs=wcs)
    >>> plt.imshow(cutout.data, origin='lower')

.. plot::

    import numpy as np
    from astropy.modeling.models import Gaussian2D
    from astropy.nddata import Cutout
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
    wcs.wcs.crpix = [100, 50]
    cutout = Cutout(data, position, (30, 40), wcs=wcs)
    plt.imshow(cutout.data)

The ``wcs`` attribute of the `~astropy.nddata.utils.Cutout` object now
contains the propagated `~astropy.wcs.WCS` for the cutout array.
Let's find the sky coordinates for a given pixel in the cutout array.
Note that we need to use the ``cutout.wcs`` object for the cutout
positions::

    >>> from astropy.wcs.utils import pixel_to_skycoord
    >>> x_cutout, y_cutout = (5, 10)
    >>> pixel_to_skycoord(x_cutout, y_cutout, cutout.wcs)
    <SkyCoord (ICRS): (ra, dec) in deg
        (197.87384041, -1.32233044)>

We now find the corresponding pixel in the large ``data`` array and
its sky coordinates::

    >>> y_data, x_data = cutout.to_large((y_cutout, x_cutout))
    >>> pixel_to_skycoord(x_data, y_data, wcs)
    <SkyCoord (ICRS): (ra, dec) in deg
        (197.87384041, -1.32233044)>

As expected, the sky coordinates in the original ``data`` and the
cutout array agree.


Reference/API
=============

.. automodapi:: astropy.nddata.utils
    :no-inheritance-diagram:
