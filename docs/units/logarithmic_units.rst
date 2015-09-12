.. |quantity| replace:: :class:`~astropy.units.Quantity`

.. _logarithmic_units:

Magnitudes and other Logarithmic Units
======================================

Magnitudes and logarithmic units such as ``dex`` and ``dB`` are used the
logarithm of values relative to some reference value.  Quantities with such
units are supported in ``astropy`` via the :class:`~astropy.units.Magnitude`,
:class:`~astropy.units.Dex`, and :class:`~astropy.units.Decibel` classes.

Creating Logarithmic Quantities
-------------------------------

One can create logarithmic quantities either directly or by multiplication with
a logarithmic unit.  For instance::

  >>> import astropy.units as u, astropy.constants as c, numpy as np
  >>> u.Magnitude(-10.)
  <Magnitude -10.0 mag>
  >>> u.Magnitude(10 * u.ct / u.s)
  <Magnitude -2.5 mag(ct / s)>
  >>> u.Magnitude(-2.5, "mag(ct/s)")
  <Magnitude -2.5 mag(ct / s)>
  >>> -2.5 * u.mag(u.ct / u.s)
  <Magnitude -2.5 mag(ct / s)>
  >>> u.Dex((c.G * u.M_sun / u.R_sun**2).cgs)  # doctest: +FLOAT_CMP
  <Dex 4.438443765928838 dex(cm / s2)>
  >>> np.linspace(2., 5., 7) * u.Unit("dex(cm/s2)")
  <Dex [ 2. , 2.5, 3. , 3.5, 4. , 4.5, 5. ] dex(cm / s2)>

Above, we make use of the fact that the units ``mag``, ``dex``, and
``dB`` are special in that, when used as functions, they return a
:class:`~astropy.units.function.logarithmic.LogUnit` instance
(:class:`~astropy.units.function.logarithmic.MagUnit`,
:class:`~astropy.units.function.logarithmic.DexUnit`, and
:class:`~astropy.units.function.logarithmic.DecibelUnit`,
respectively).  The same happens as required when strings are parsed
by :class:`~astropy.units.Unit`.

As for normal |quantity| objects, one can access the value with the
`~astropy.units.Quantity.value` attribute. In addition, one can convert easily
to a |quantity| with the physical unit using the
`~astropy.units.function.FunctionQuantity.physical` attribute::

    >>> logg = 5. * u.dex(u.cm / u.s**2)
    >>> logg.value
    5.0
    >>> logg.physical
    <Quantity 100000.0 cm / s2>

Converting to different units
-----------------------------

Like |quantity| objects, logarithmic quantities can be converted to different
units using the :meth:`~astropy.units.function.FunctionQuantity.to` method.
Here, if the requested unit is not a logarithmic unit, the object will be
automatically converted to its physical unit::

    >>> logg = 5. * u.dex(u.cm / u.s**2)
    >>> logg.to(u.m / u.s**2)
    <Quantity 1000.0 m / s2>
    >>> logg.to('dex(m/s2)')
    <Dex 3.0 dex(m / s2)>

For convenience, the `~astropy.units.function.FunctionQuantity.si` and
`~astropy.units.function.FunctionQuantity.cgs` attributes can be used
to convert the |quantity| to base S.I. or c.g.s units::

    >>> logg.si
    <Dex 3.0 dex(m / s2)>

Arithmetic
----------

Addition and subtraction work as expected for logarithmic quantities,
multiplying and dividing the physical units as appropriate.  It may be best
seen through an example of a very simple photometric reduction.  First,
calculate instrumental magnitudes assuming some count rates for three objects::

    >>> tint = 1000.*u.s
    >>> cr_b = ([3000., 100., 15.] * u.ct) / tint
    >>> cr_v = ([4000., 90., 25.] * u.ct) / tint
    >>> b_i, v_i = u.Magnitude(cr_b), u.Magnitude(cr_v)
    >>> b_i, v_i  # doctest: +FLOAT_CMP
    (<Magnitude [-1.19280314, 2.5       , 4.55977185] mag(ct / s)>,
     <Magnitude [-1.50514998, 2.61439373, 4.00514998] mag(ct / s)>)

Then, the instrumental B-V color is simply::

    >>> b_i - v_i
    <Magnitude [ 0.31234684,-0.11439373, 0.55462187] mag>

Note that the physical unit has become dimensionless.  The following step might
be used to correct for atmospheric extinction::

    >>> atm_ext_b, atm_ext_v = 0.12 * u.mag, 0.08 * u.mag
    >>> secz = 1./np.cos(45 * u.deg)
    >>> b_i0 = b_i - atm_ext_b * secz
    >>> v_i0 = v_i - atm_ext_b * secz
    >>> b_i0, v_i0  # doctest: +FLOAT_CMP
    (<Magnitude [-1.36250876, 2.33029437, 4.39006622] mag(ct / s)>,
     <Magnitude [-1.67485561, 2.4446881 , 3.83544435] mag(ct / s)>)

Since the extinction is dimensionless, the units do not change.  Now suppose
the first star has a known ST magnitude, so we can calculate zero points::

    >>> b_ref, v_ref = 17.2 * u.STmag, 17.0 * u.STmag
    >>> b_ref, v_ref  # doctest: +FLOAT_CMP
    (<Magnitude 17.2 mag(ST)>, <Magnitude 17.0 mag(ST)>)
    >>> zp_b, zp_v = b_ref - b_i0[0], v_ref - v_i0[0]
    >>> zp_b, zp_v  # doctest: +FLOAT_CMP
    (<Magnitude 18.562508764283926 mag(s ST / ct)>,
     <Magnitude 18.674855605804677 mag(s ST / ct)>)

Here, ``ST`` is a short-hand for the ST zero-point flux::

    >>> (0. * u.STmag).to(u.erg/u.s/u.cm**2/u.AA)  # doctest: +FLOAT_CMP
    <Quantity 3.6307805477010028e-09 erg / (Angstrom cm2 s)>
    >>> (-21.1 * u.STmag).to(u.erg/u.s/u.cm**2/u.AA)  # doctest: +FLOAT_CMP
    <Quantity 1. erg / (Angstrom cm2 s)>

.. note:: only ST [H+95]_ and AB [OG83]_ magnitudes are implemented at
	  present, as these are defined in terms of flux densities, i.e.,
          do not depend on the filter the measurement was made with.

Now applying the calibration, we find (note the proper change in units)::

    >>> B, V = b_i0 + zp_b, v_i0 + zp_v
    >>> B, V  # doctest: +FLOAT_CMP
    (<Magnitude [ 17.2       , 20.89280314, 22.95257499] mag(ST)>,
     <Magnitude [ 17.        , 21.1195437 , 22.51029996] mag(ST)>)

We could convert these magnitudes to another system, e.g., ABMag, using
appropriate equivalency::

    >>> V.to(u.ABmag, u.spectral_density(5500.*u.AA))  # doctest: +FLOAT_CMP
    <Magnitude [ 16.99023831, 21.10978201, 22.50053827] mag(AB)>

Suppose we also knew the intrinsic color of the first star, then we can
calculate the reddening::

    >>> B_V0 = -0.2 * u.mag
    >>> EB_V = (B - V)[0] - B_V0
    >>> R_V = 3.1
    >>> A_V = R_V * EB_V
    >>> A_B = (R_V+1) * EB_V
    >>> EB_V, A_V, A_B  # doctest: +FLOAT_CMP
    (<Magnitude 0.3999999999999993 mag>,
     <Quantity 1.2399999999999978 mag>,
     <Quantity 1.639999999999997 mag>)

Here, one sees that the extinctions have been converted to quantities. This
happens generally for division and multiplication, since these processes
work only for dimensionless magnitudes (otherwise, the physical unit would have
to be raised to some power), and |quantity| objects, unlike logarithmic
quantities, allow units like ``mag / d``.

Note that one can take the automatic unit conversion quite far (perhaps too
far, but it is fun).  For instance, suppose we also knew the absolute
magnitude, then we can define the appropriate corresponding luminosity and
absolute magnitude and calculate the distance modulus::

    >>> ST0abs = u.Unit('STabs', u.STmag.physical_unit * 4.*np.pi*(10.*u.pc)**2)
    >>> STabsmag = u.mag(ST0abs)
    >>> M_V = 5.76 * STabsmag
    >>> M_B = M_V + B_V0
    >>> DM = V[0] - A_V - M_V
    >>> M_V, M_B, DM  # doctest: +FLOAT_CMP
    (<Magnitude 5.76 mag(STabs)>,
     <Magnitude 5.56 mag(STabs)>,
     <Magnitude 10.000000000000002 mag(ST / STabs)>)

With a proper equivalency, we can also convert to distance without remembering
the 5-5log rule::

    >>> radius_and_inverse_area = [(u.pc, u.pc**-2,
    ...                            lambda x: 1./(4.*np.pi*x**2),
    ...                            lambda x: np.sqrt(1./(4.*np.pi*x)))]
    >>> DM.to(u.pc, equivalencies=radius_and_inverse_area)  # doctest: +FLOAT_CMP
    <Quantity 1000.0000000000009 pc>

Numpy functions
---------------

For logarithmic quantities, most numpy functions and many array methods do not
make sense, hence they are disabled.  But one can use those one would expect to
work::

    >>> np.max(v_i)  # doctest: +FLOAT_CMP
    <Magnitude 4.005149978319905 mag(ct / s)>
    >>> np.std(v_i)  # doctest: +FLOAT_CMP
    <Magnitude 2.339711494548601 mag>

.. note:: This is implemented by having a list of supported ufuncs in
	  ``units/function/core.py`` and by explicitly disabling some
	  array methods in :class:`~astropy.units.function.FunctionQuantity`.
          If you believe a function or method is incorrectly treated,
	  please `let us know <http://www.astropy.org/contribute.html>`_.

Dimensionless logarithmic quantities
------------------------------------

Dimensionless quantities are treated somewhat specially, in that, if needed,
logarithmic quantities will be converted to normal |quantity| objects with the
appropriate unit of ``mag``, ``dB``, or ``dex``.  With this, it is possible to
use composite units like ``mag/d`` or ``dB/m``, which cannot easily be
supported as logarithmic units.  For instance::

    >>> dBm = u.dB(u.mW)
    >>> signal_in, signal_out = 100. * dBm, 50 * dBm
    >>> cable_loss = (signal_in - signal_out) / (100. * u.m)
    >>> signal_in, signal_out, cable_loss
    (<Decibel 100.0 dB(mW)>, <Decibel 50.0 dB(mW)>, <Quantity 0.5 dB / m>)
    >>> better_cable_loss = 0.2 * u.dB / u.m
    >>> signal_in - better_cable_loss * 100. * u.m
    <Decibel 80.0 dB(mW)>



.. [H+95] E.g., Holtzman et al., 1995, `PASP 107, 1065
          <http://adsabs.harvard.edu/abs/1995PASP..107.1065H>`_
.. [OG83] Oke, J.B., & Gunn, J. E., 1983, `ApJ 266, 713
	  <http://adsabs.harvard.edu/abs/1983ApJ...266..713O>`_
