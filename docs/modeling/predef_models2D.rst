.. _predef_models2D:

*********
2D Models
*********

These models take as input x and y arrays.

Operations
==========

These models perform simple mathematical operations.

- :class:`~astropy.modeling.functional_models.Const2D` model returns the
  constant replicated by the number of input x and y values.

Shapes
======

These models provide shapes, often used to model general x, y, z data.

- :class:`~astropy.modeling.functional_models.Planar2D` model computes
  a tilted plan with specified x,y slopes and z intercept

Profiles
========

These models provide profiles, often used sources in images.
All models have parameters giving the x,y location of the center and
an amplitude.

- :class:`~astropy.modeling.functional_models.AiryDisk2D` model computes
  the Airy function for a radius

- :class:`~astropy.modeling.functional_models.Box2D` model computes a box
  with x,y dimensions

- :class:`~astropy.modeling.functional_models.Disk2D` model computes a
  disk a radius

- :class:`~astropy.modeling.functional_models.Ellipse2D` model computes
  an ellipse with major and minor axis and rotation angle

- :class:`~astropy.modeling.functional_models.Gaussian2D` model computes
  a Gaussian with x,y standard deviations and rotation angle

- :class:`~astropy.modeling.functional_models.RickerWavelet2D` model computes
  a symmetric RickerWavelet function with the specified sigma

- :class:`~astropy.modeling.functional_models.Sersic2D` model computes
  a Sersic profile with an effective half-light radius, rotation, and
  Sersic index

- :class:`~astropy.modeling.functional_models.TrapezoidDisk2D` model
  computes a disk with a radius and slope

- :class:`~astropy.modeling.functional_models.Ring2D` model computes
  a ring with inner and outer radii

.. plot::

    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    from astropy.modeling.models import (AiryDisk2D, Box2D, Disk2D, Ellipse2D,
                                         Gaussian2D, RickerWavelet2D, Sersic2D,
                                         TrapezoidDisk2D, Ring2D)

    x = np.linspace(-4.0, 6.0, num=100)
    r = np.logspace(-1.0, 2.0, num=100)

    fig, sax = plt.subplots(nrows=3, ncols=3, figsize=(8, 8))
    ax = sax.flatten()

    # setup the x,y coordinates
    x_npts = 100
    y_npts = x_npts
    x0, x1 = -4, 6
    y0, y1 = -3, 7
    x = np.linspace(x0, x1, num=x_npts)
    y = np.linspace(y0, y1, num=y_npts)
    X, Y = np.meshgrid(x, y)

    # plot the different 2D profiles
    mods = [AiryDisk2D(amplitude=10.0, x_0=1.0, y_0=2.0, radius=1.0),
            Box2D(amplitude=10.0, x_0=1.0, y_0=2.0, x_width=1.0, y_width=2.0),
            Disk2D(amplitude=10.0, x_0=1.0, y_0=2.0, R_0=1.0),
            Ellipse2D(amplitude=10.0, x_0=1.0, y_0=2.0, a=1.0, b=2.0, theta=math.pi/4.),
            Gaussian2D(amplitude=10.0, x_mean=1.0, y_mean=2.0, x_stddev=1.0, y_stddev=2.0, theta=math.pi/4.),
            RickerWavelet2D(amplitude=10.0, x_0=1.0, y_0=2.0, sigma=1.0),
            Sersic2D(amplitude=10.0, x_0=1.0, y_0=2.0, r_eff=1.0, ellip=0.5, theta=math.pi/4.),
            TrapezoidDisk2D(amplitude=10.0, x_0=1.0, y_0=2.0, R_0=1.0, slope=5.0),
            Ring2D(amplitude=10.0, x_0=1.0, y_0=2.0, r_in=1.0, r_out=2.0)]

    for k, mod in enumerate(mods):
        cname = mod.__class__.__name__
        ax[k].set_title(cname)
        if cname == "AiryDisk2D":
            normfunc = LogNorm(vmin=0.001, vmax=10.)
        elif cname in ["Gaussian2D", "Sersic2D"]:
            normfunc = LogNorm(vmin=0.1, vmax=10.)
        else:
            normfunc = None
        ax[k].imshow(mod(X, Y), extent=[x0, x1, y0, y1], origin="lower", cmap=plt.cm.gray_r,
                     norm=normfunc)

    for k in range(len(mods)):
        ax[k].set_xlabel("x")
        ax[k].set_ylabel("y")

    # remove axis for any plots not used
    for k in range(len(mods), len(ax)):
        ax[k].axis("off")

    plt.tight_layout()
    plt.show()
