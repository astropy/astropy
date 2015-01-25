.. _coordinates-galactocentric:

==========================================
Transforming to Galactocentric coordinates
==========================================

This document describes the mathematics behind the transformation from
:class:`~astropy.coordinates.ICRS` to `~astropy.coordinates.Galactocentric`
coordinates. For examples of how to use this transformation in code, see
the *Examples* section of the
`~astropy.coordinates.Galactocentric` class documentation.

We assume that we start with a 3D position in the ICRS reference frame:
a Right Ascension, Declination, and heliocentric distance,
:math:`(\alpha, \delta, d)`. We can trivially convert this to a
Cartesian position using the standard transformation from Cartesian to
spherical coordinates:

.. math::

   \begin{aligned}
       x_{\rm icrs} &= d\cos{\alpha}\cos{\delta}\\
       y_{\rm icrs} &= d\sin{\alpha}\cos{\delta}\\
       z_{\rm icrs} &= d\sin{\delta}\\
       \boldsymbol{r}_{\rm icrs} &= \begin{pmatrix}
         x_{\rm icrs}\\
         y_{\rm icrs}\\
         z_{\rm icrs}
       \end{pmatrix}\end{aligned}

The first transformations will rotate the :math:`x_{\rm icrs}` axis so
that the new :math:`x'` axis points towards the Galactic Center (GC),
specified by the ICRS position
:math:`(\alpha_{\rm GC}, \delta_{\rm GC})`:

.. math::

   \begin{aligned}
       \boldsymbol{R}_1 &= \begin{bmatrix}
         \cos\delta_{\rm GC}& 0 & -\sin\delta_{\rm GC}\\
         0 & 1 & 0 \\
         \sin\delta_{\rm GC}& 0 & \cos\delta_{\rm GC}\end{bmatrix}\\
       \boldsymbol{R}_2 &=
       \begin{bmatrix}
         \cos\alpha_{\rm GC}& \sin\alpha_{\rm GC}& 0\\
         -\sin\alpha_{\rm GC}& \cos\alpha_{\rm GC}& 0\\
         0 & 0 & 1
       \end{bmatrix}.\end{aligned}

The transformation thus far has aligned the :math:`x'` axis with the
vector pointing from the Sun to the GC, but the :math:`y'` and
:math:`z'` axes point in an arbitrary direction. We adopt the
orientation of the Galactic plane as the normal to the north pole of
Galactic coordinates defined by the IAU
(`Blaauw et. al. 1960 <http://adsabs.harvard.edu/abs/1960MNRAS.121..164B>`_).
This extra “roll” angle, :math:`\eta`, was measured by transforming a grid
of points along :math:`l=0` to this interim frame and minimizing the square
of their :math:`y'` positions. We find:

.. math::

   \begin{aligned}
       \eta &= 148.5986320^\circ\\
       \boldsymbol{R}_3 &=
       \begin{bmatrix}
         1 & 0 & 0\\
         0 & \cos\eta & \sin\eta\\
         0 & -\sin\eta & \cos\eta
       \end{bmatrix}\end{aligned}

The full rotation matrix thus far is:

.. math::

   \begin{gathered}
       \boldsymbol{R} = \boldsymbol{R}_3 \boldsymbol{R}_1 \boldsymbol{R}_2 = \\
       \begin{bmatrix}
         \cos\alpha_{\rm GC}\cos\delta_{\rm GC}& \cos\delta_{\rm GC}\sin\alpha_{\rm GC}& -\sin\delta_{\rm GC}\\
         \cos\alpha_{\rm GC}\sin\delta_{\rm GC}\sin\eta - \sin\alpha_{\rm GC}\cos\eta & \sin\alpha_{\rm GC}\sin\delta_{\rm GC}\sin\eta + \cos\alpha_{\rm GC}\cos\eta & \cos\delta_{\rm GC}\sin\eta\\
         \cos\alpha_{\rm GC}\sin\delta_{\rm GC}\cos\eta + \sin\alpha_{\rm GC}\sin\eta & \sin\alpha_{\rm GC}\sin\delta_{\rm GC}\cos\eta - \cos\alpha_{\rm GC}\sin\eta & \cos\delta_{\rm GC}\cos\eta
       \end{bmatrix}\end{gathered}

With the rotated position vector
:math:`\boldsymbol{R}\boldsymbol{r}_{\rm icrs}`, we can now subtract the
distance to the GC, :math:`d_{\rm GC}`, which is purely along the
:math:`x'` axis:

.. math::

   \begin{aligned}
       \boldsymbol{r}' &= \boldsymbol{R}\boldsymbol{r}_{\rm icrs} - d_{\rm GC}\hat{\boldsymbol{x}}_{\rm GC}.\end{aligned}

The final transformation is to account for the height of the Sun above
the Galactic midplane by rotating about the final :math:`y''` axis by
the angle :math:`\theta= \sin^{-1}(z_\odot / d_{\rm GC})`:

.. math::

   \begin{aligned}
       \boldsymbol{H} &=
       \begin{bmatrix}
         \cos\theta & 0 & \sin\theta\\
         0 & 1 & 0\\
         -\sin\theta & 0 & \cos\theta
       \end{bmatrix}\end{aligned}

where :math:`z_\odot` is the measured height of the Sun above the
midplane.

The full transformation is then:

.. math:: \boldsymbol{r}_{\rm GC} = \boldsymbol{H} \left( \boldsymbol{R}\boldsymbol{r}_{\rm icrs} - d_{\rm GC}\hat{\boldsymbol{x}}_{\rm GC}\right).

