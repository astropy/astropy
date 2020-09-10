.. _coordinates-galactocentric:

**************************************************
Description of the Galactocentric Coordinate Frame
**************************************************

While many other frames implemented in `astropy.coordinates` are standardized in
some way (e.g., defined by the IAU), there is no standard Milky Way
reference frame with the center of the Milky Way as its origin. (This is
distinct from `~astropy.coordinates.Galactic` coordinates, which point
toward the Galactic Center but have their origin in the Solar System).
The `~astropy.coordinates.Galactocentric` frame
class is meant to be flexible enough to support all common definitions of such a
transformation, but with reasonable default parameter values, such as the solar
velocity relative to the Galactic center, the solar height above the Galactic
midplane, etc. Below, `we describe our generalized definition of the
transformation <astropy-coordinates-galactocentric-transformation>`_ from the
ICRS to/from Galactocentric coordinates, and `describe how to customize the
default Galactocentric parameters
<astropy-coordinates-galactocentric-defaults>`_ that are used when the
`~astropy.coordinates.Galactocentric` frame is initialized without explicitly
passing in parameter values.


.. _astropy-coordinates-galactocentric-transformation:

Definition of the Transformation
================================

This document describes the mathematics behind the transformation from
`~astropy.coordinates.ICRS` to `~astropy.coordinates.Galactocentric`
coordinates. This is described in detail here on account of the mathematical
subtleties and the fact that there is no official standard/definition for this
frame. For examples of how to use this transformation in code, see the
the *Examples* section of the `~astropy.coordinates.Galactocentric` class
documentation.

We assume that we start with a 3D position in the ICRS reference frame:
a Right Ascension, Declination, and heliocentric distance,
:math:`(\alpha, \delta, d)`. We can convert this to a Cartesian position using
the standard transformation from Cartesian to spherical coordinates:

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

The first transformation rotates the :math:`x_{\rm icrs}` axis so that the new
:math:`x'` axis points towards the Galactic Center (GC), specified by the ICRS
position :math:`(\alpha_{\rm GC}, \delta_{\rm GC})` (in the
`~astropy.coordinates.Galactocentric` frame, this is controlled by the frame
attribute ``galcen_coord``):

.. math::

   \begin{aligned}
       \boldsymbol{R}_1 &= \begin{bmatrix}
         \cos\delta_{\rm GC}& 0 & \sin\delta_{\rm GC}\\
         0 & 1 & 0 \\
         -\sin\delta_{\rm GC}& 0 & \cos\delta_{\rm GC}\end{bmatrix}\\
       \boldsymbol{R}_2 &=
       \begin{bmatrix}
         \cos\alpha_{\rm GC}& \sin\alpha_{\rm GC}& 0\\
         -\sin\alpha_{\rm GC}& \cos\alpha_{\rm GC}& 0\\
         0 & 0 & 1
       \end{bmatrix}.\end{aligned}

The transformation thus far has aligned the :math:`x'` axis with the
vector pointing from the Sun to the GC, but the :math:`y'` and
:math:`z'` axes point in arbitrary directions. We adopt the
orientation of the Galactic plane as the normal to the north pole of
Galactic coordinates defined by the IAU
(`Blaauw et. al. 1960 <https://ui.adsabs.harvard.edu/abs/1960MNRAS.121..164B>`_).
This extra “roll” angle, :math:`\eta`, was measured by transforming a grid
of points along :math:`l=0` to this interim frame and minimizing the square
of their :math:`y'` positions. We find:

.. math::

   \begin{aligned}
       \eta &= 58.5986320306^\circ\\
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

where :math:`\hat{\boldsymbol{x}}_{\rm GC} = (1,0,0)^{\mathsf{T}}`.

The final transformation accounts for the (specified) height of the Sun above
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

.. topic:: Examples:

    For an example of how to use the `~astropy.coordinates.Galactocentric`
    frame, see
    :ref:`sphx_glr_generated_examples_coordinates_plot_galactocentric-frame.py`.


.. _astropy-coordinates-galactocentric-defaults:

Controlling the Default Frame Parameters
========================================

All of the frame-defining parameters of the
`~astropy.coordinates.Galactocentric` frame are customizable and can be set by
passing arguments in to the `~astropy.coordinates.Galactocentric` initializer.
However, it is often convenient to use the frame without having to pass in every
parameter. Hence, the class comes with reasonable default values for these
parameters, but more precise measurements of the solar position or motion in the
Galaxy are constantly being made. The default values of the
`~astropy.coordinates.Galactocentric` frame attributes will therefore be updated
as necessary with subsequent releases of ``astropy``. We therefore provide a
mechanism to globally or locally control the default parameter values used in
this frame through the `~astropy.coordinates.galactocentric_frame_defaults`
`~astropy.utils.state.ScienceState` class.

The `~astropy.coordinates.galactocentric_frame_defaults` class controls the
default parameter settings in `~astropy.coordinates.Galactocentric` by mapping a
set of string names to particular choices of the parameter values. For an
up-to-date list of valid names, see the docstring of
`~astropy.coordinates.galactocentric_frame_defaults`, but these names are things
like ``'pre-v4.0'``, which sets the default parameter values to their original
definition (i.e. pre-astropy-v4.0) values, and ``'v4.0'``, which sets the
default parameter values to a more modern set of measurements as updated in
Astropy version 4.0.

As with other `~astropy.utils.state.ScienceState` subclasses, the
`~astropy.coordinates.galactocentric_frame_defaults` class can be used to
globally set the frame defaults at runtime. For example, the default parameter
values can be seen by initializing the `~astropy.coordinates.Galactocentric`
frame with no arguments:

.. testsetup::

    >>> from astropy.coordinates import galactocentric_frame_defaults
    >>> _ = galactocentric_frame_defaults.set('pre-v4.0')

::

    >>> from astropy.coordinates import Galactocentric
    >>> Galactocentric()
    <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
        (266.4051, -28.936175)>, galcen_distance=8.3 kpc, galcen_v_sun=(11.1, 232.24, 7.25) km / s, z_sun=27.0 pc, roll=0.0 deg)>

These default values can be modified using this class::

    >>> from astropy.coordinates import galactocentric_frame_defaults
    >>> _ = galactocentric_frame_defaults.set('v4.0')
    >>> Galactocentric() # doctest: +FLOAT_CMP
    <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
        (266.4051, -28.936175)>, galcen_distance=8.122 kpc, galcen_v_sun=(12.9, 245.6, 7.78) km / s, z_sun=20.8 pc, roll=0.0 deg)>
    >>> _ = galactocentric_frame_defaults.set('pre-v4.0')
    >>> Galactocentric() # doctest: +FLOAT_CMP
    <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
        (266.4051, -28.936175)>, galcen_distance=8.3 kpc, galcen_v_sun=(11.1, 232.24, 7.25) km / s, z_sun=27.0 pc, roll=0.0 deg)>

The default parameters can also be updated by using this class as a context
manager to change the default parameter values locally to a piece of your code::

    >>> with galactocentric_frame_defaults.set('pre-v4.0'):
    ...     print(Galactocentric()) # doctest: +FLOAT_CMP
    <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
        (266.4051, -28.936175)>, galcen_distance=8.3 kpc, galcen_v_sun=(11.1, 232.24, 7.25) km / s, z_sun=27.0 pc, roll=0.0 deg)>

Again, changing the default parameter values will not affect frame
attributes that are explicitly specified::

    >>> import astropy.units as u
    >>> with galactocentric_frame_defaults.set('pre-v4.0'):
    ...     print(Galactocentric(galcen_distance=8.0*u.kpc)) # doctest: +FLOAT_CMP
    <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
        (266.4051, -28.936175)>, galcen_distance=8.0 kpc, galcen_v_sun=(11.1, 232.24, 7.25) km / s, z_sun=27.0 pc, roll=0.0 deg)>

Starting with Astropy v4.1, unless set with the
`~astropy.coordinates.galactocentric_frame_defaults` class, the default
parameter values for the `~astropy.coordinates.Galactocentric` frame will be set
to ``'latest'``, meaning that the default parameter values may change if you
update Astropy. If you use the `~astropy.coordinates.Galactocentric` frame
without specifying all parameter values explicitly, we therefore suggest
manually setting the frame default set manually in any science code that depends
sensitively on the choice of, e.g., solar motion or the other frame parameters.
For example, in such code, we recommend adding something like this to your
import block (here using ``'v4.0'`` as an example)::

    >>> import astropy.coordinates as coord
    >>> coord.galactocentric_frame_defaults.set('v4.0') # doctest: +SKIP
