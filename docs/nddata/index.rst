.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

The `~astropy.nddata` package provides classes to represent images and other gridded data, some essential functions for manipulating images, and the infrastructure for package developers who wish to include support for the image classes.

Getting started
===============

Though the `~astropy.nddata` package supports any kind of gridded data, this introduction will focus on the use of `~astropy.nddata` for two-dimensional images. To get started, we'll construct a two-dimensional image with a few sources, some Gaussian noise, and a cosmic ray which we will later mask out::

    >>> import numpy as np
    >>> from astropy.modeling.models import Gaussian2D
    >>> y, x = np.mgrid[0:500, 0:600]
    >>> data = (Gaussian2D(1, 50, 100, 20, 10, theta=0.5)(x, y) +
    ...         Gaussian2D(0.5, 400, 300, 8, 12, theta=1.2)(x,y) +
    ...         Gaussian2D(0.75, 250, 400, 5, 7, theta=0.23)(x,y) +
    ...         Gaussian2D(0.9, 525, 150, 3, 3)(x,y) +
    ...         Gaussian2D(0.6, 200, 225, 3, 3)(x,y))
    >>> data += 0.01 * np.random.randn(500, 600)
    >>> cosmic_ray_value = 0.997
    >>> data[300:310, 100] = cosmic_ray_value

This image has a large "galaxy" in the lower left and the cosmic ray is the horizontal line in the lower middle of the image:

.. doctest-skip::

    >>> import matplotlib.pyplot as plt
    >>> plt.imshow(data, origin='lower')

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    y, x = np.mgrid[0:500, 0:600]
    data = (Gaussian2D(1, 50, 100, 20, 10, theta=0.5)(x, y) +
            Gaussian2D(0.5, 400, 300, 8, 12, theta=1.2)(x,y) +
            Gaussian2D(0.75, 250, 400, 5, 7, theta=0.23)(x,y) +
            Gaussian2D(0.9, 525, 150, 3, 3)(x,y) +
            Gaussian2D(0.6, 200, 225, 3, 3)(x,y))
    data += 0.01 * np.random.randn(500, 600)
    cosmic_ray_value = 0.997
    data[300:310, 100] = cosmic_ray_value
    plt.imshow(data, origin='lower')


Foo!




Using ``nddata``
================

.. toctree::
   :maxdepth: 2

   nddata.rst
   decorator.rst
   mixins/index.rst
   subclassing.rst
   utils.rst

Reference/API
=============


.. _APE 7: https://github.com/astropy/astropy-APEs/blob/master/APE7.rst
