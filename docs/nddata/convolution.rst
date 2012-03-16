Convolution
===========

Introduction
------------

``astropy.nddata`` includes a convolution function that offers improvements compared to the scipy ``astropy.ndimage`` convolution routines, including:

* Proper treatment of NaN values

* A single function for 1-D, 2-D, and 3-D convolution

* Improved options for the treatment of edges

The following thumbnails show the difference between Scipy's and Astropy's convolve functions on an Astronomical image that contains NaN values. Scipy's function essentially returns NaN for all pixels that are within a kernel of any NaN value, which is often not the desired result.

.. |original| image:: images/original.png
.. |scipy| image:: images/scipy.png
.. |astropy| image:: images/astropy.png

+-----------------------+--------------------+----------------------+
|        Original       | Scipy ``convolve`` | Astropy ``convolve`` |
+-----------------------+--------------------+----------------------+
|       |original|      |       |scipy|      |      |astropy|       | 
+-----------------------+--------------------+----------------------+


Usage
-----

The convolution function is imported with::

    from astropy.nddata import convolve
    
and is used as::

    result = convolve(image, kernel)

The input images and kernels should be lists or Numpy arrays with either both 1, 2, or 3 dimensions (and the number of dimensions should be the same for the image and kernel). The result is a Numpy array with the same dimensions as the input image.

The ``convolve`` function takes an optional ``boundary=`` argument describing how to perform the convolution at the edge of the array. The values for ``boundary`` can be:

* ``None``: set the result values to zero where the kernel extends beyond the edge of the array (default)

* ``'fill'``: set values outside the array boundary to a constant. If this option is specified, the constant should be specified using the ``fill_value=`` argument, which defaults to zero.

* ``'wrap'``: assume that the boundaries are periodic

* ``'extend'`` : set values outside the array to the nearest array value

By default, the kernel is not normalized. To normalize it prior to convolution, use::

    result = convolve(image, kernel, normalize_kernel=True)