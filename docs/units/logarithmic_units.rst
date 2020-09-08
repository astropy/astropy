.. |quantity| replace:: :class:`~astropy.units.Quantity`

.. _logarithmic_units:

Magnitudes and other Logarithmic Units
**************************************

Magnitudes and logarithmic units such as ``dex`` and ``dB`` are used the
logarithm of values relative to some reference value.  Quantities with such
units are supported in ``astropy`` via the :class:`~astropy.units.Magnitude`,
:class:`~astropy.units.Dex`, and :class:`~astropy.units.Decibel` classes.

Creating Logarithmic Quantities
===============================

One can create logarithmic quantities either directly or by multiplication with
a logarithmic unit.  For instance::

  >>> import astropy.units as u, astropy.constants as c, numpy as np
  >>> u.Magnitude(-10.)  # doctest: +FLOAT_CMP
  <Magnitude -10. mag>
  >>> u.Magnitude(10 * u.ct / u.s)  # doctest: +FLOAT_CMP
  <Magnitude -2.5 mag(ct / s)>
  >>> u.Magnitude(-2.5, "mag(ct/s)")  # doctest: +FLOAT_CMP
  <Magnitude -2.5 mag(ct / s)>
  >>> -2.5 * u.mag(u.ct / u.s)  # doctest: +FLOAT_CMP
  <Magnitude -2.5 mag(ct / s)>
  >>> u.Dex((c.G * u.M_sun / u.R_sun**2).cgs)  # doctest: +FLOAT_CMP
  <Dex 4.438067627303133 dex(cm / s2)>
  >>> np.linspace(2., 5., 7) * u.Unit("dex(cm/s2)")  # doctest: +FLOAT_CMP
  <Dex [2. , 2.5, 3. , 3.5, 4. , 4.5, 5. ] dex(cm / s2)>

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
    >>> logg.physical  # doctest: +FLOAT_CMP
    <Quantity 100000. cm / s2>

Converting to different units
=============================

Like |quantity| objects, logarithmic quantities can be converted to different
units, be it another logarithmic unit or a physical one::

    >>> logg = 5. * u.dex(u.cm / u.s**2)
    >>> logg.to(u.m / u.s**2)  # doctest: +FLOAT_CMP
    <Quantity 1000. m / s2>
    >>> logg.to('dex(m/s2)')  # doctest: +FLOAT_CMP
    <Dex 3. dex(m / s2)>

For convenience, the `~astropy.units.function.FunctionQuantity.si` and
`~astropy.units.function.FunctionQuantity.cgs` attributes can be used
to convert the |quantity| to base S.I. or c.g.s units::

    >>> logg.si  # doctest: +FLOAT_CMP
    <Dex 3. dex(m / s2)>

Arithmetic and Photometric Applications
=======================================

Addition and subtraction work as expected for logarithmic quantities,
multiplying and dividing the physical units as appropriate.  It may be best
seen through an example of a very simple photometric reduction.  First,
calculate instrumental magnitudes assuming some count rates for three objects::

    >>> tint = 1000.*u.s
    >>> cr_b = ([3000., 100., 15.] * u.ct) / tint
    >>> cr_v = ([4000., 90., 25.] * u.ct) / tint
    >>> b_i, v_i = u.Magnitude(cr_b), u.Magnitude(cr_v)
    >>> b_i, v_i  # doctest: +FLOAT_CMP
    (<Magnitude [-1.19280314,  2.5       ,  4.55977185] mag(ct / s)>,
     <Magnitude [-1.50514998,  2.61439373,  4.00514998] mag(ct / s)>)

Then, the instrumental B-V color is simply::

    >>> b_i - v_i  # doctest: +FLOAT_CMP
    <Magnitude [ 0.31234684, -0.11439373,  0.55462187] mag>

Note that the physical unit has become dimensionless.  The following step might
be used to correct for atmospheric extinction::

    >>> atm_ext_b, atm_ext_v = 0.12 * u.mag, 0.08 * u.mag
    >>> secz = 1./np.cos(45 * u.deg)
    >>> b_i0 = b_i - atm_ext_b * secz
    >>> v_i0 = v_i - atm_ext_b * secz
    >>> b_i0, v_i0  # doctest: +FLOAT_CMP
    (<Magnitude [-1.36250876,  2.33029437,  4.39006622] mag(ct / s)>,
     <Magnitude [-1.67485561,  2.4446881 ,  3.83544435] mag(ct / s)>)

Since the extinction is dimensionless, the units do not change.  Now suppose
the first star has a known ST magnitude, so we can calculate zero points::

    >>> b_ref, v_ref = 17.2 * u.STmag, 17.0 * u.STmag
    >>> b_ref, v_ref  # doctest: +FLOAT_CMP
    (<Magnitude 17.2 mag(ST)>, <Magnitude 17. mag(ST)>)
    >>> zp_b, zp_v = b_ref - b_i0[0], v_ref - v_i0[0]
    >>> zp_b, zp_v  # doctest: +FLOAT_CMP
    (<Magnitude 18.56250876 mag(s ST / ct)>,
     <Magnitude 18.67485561 mag(s ST / ct)>)

Here, ``ST`` is a short-hand for the ST zero-point flux::

    >>> (0. * u.STmag).to(u.erg/u.s/u.cm**2/u.AA)  # doctest: +FLOAT_CMP
    <Quantity 3.63078055e-09 erg / (Angstrom cm2 s)>
    >>> (-21.1 * u.STmag).to(u.erg/u.s/u.cm**2/u.AA)  # doctest: +FLOAT_CMP
    <Quantity 1. erg / (Angstrom cm2 s)>

.. note:: at present, only magnitudes defined in terms of luminosity or flux
	  are implemented, since those that do not depend on the filter the
          measurement was made with.  They include absolute and apparent
          bolometric [M15]_, ST [H95]_ and AB [OG83]_ magnitudes.

Now applying the calibration, we find (note the proper change in units)::

    >>> B, V = b_i0 + zp_b, v_i0 + zp_v
    >>> B, V  # doctest: +FLOAT_CMP
    (<Magnitude [17.2       , 20.89280314, 22.95257499] mag(ST)>,
     <Magnitude [17.        , 21.1195437 , 22.51029996] mag(ST)>)

We could convert these magnitudes to another system, e.g., ABMag, using
appropriate equivalency::

    >>> V.to(u.ABmag, u.spectral_density(5500.*u.AA))  # doctest: +FLOAT_CMP
    <Magnitude [16.99023831, 21.10978201, 22.50053827] mag(AB)>

This is particularly useful for converting magnitude into flux density.
``V`` is currently in ST magnitudes, which is based on flux densities per
unit wavelength (:math:`f_\lambda`). Therefore, we can directly convert ``V`` into 
flux density per unit wavelength using the :meth:`~astropy.units.quantity.Quantity.to` 
method::

    >>> flam = V.to(u.erg/u.s/u.cm**2/u.AA)
    >>> flam  # doctest: +FLOAT_CMP
    <Quantity [5.75439937e-16, 1.29473986e-17, 3.59649961e-18] erg / (Angstrom cm2 s)>

To convert ``V`` to flux density per unit frequency (:math:`f_\nu`), we again need 
the appropriate :ref:`equivalency <unit_equivalencies>`, which in this case is the 
central wavelength of the magnitude band, 5500 Angstroms::

    >>> lam = 5500 * u.AA
    >>> fnu = V.to(u.erg/u.s/u.cm**2/u.Hz, u.spectral_density(lam))
    >>> fnu  # doctest: +FLOAT_CMP
    <Quantity [5.80636959e-27, 1.30643316e-28, 3.62898099e-29] erg / (cm2 Hz s)>

We could have used the central frequency instead::

    >>> nu = 5.45077196e+14 * u.Hz
    >>> fnu = V.to(u.erg/u.s/u.cm**2/u.Hz, u.spectral_density(nu))
    >>> fnu  # doctest: +FLOAT_CMP
    <Quantity [5.80636959e-27, 1.30643316e-28, 3.62898099e-29] erg / (cm2 Hz s)>

.. Note::

    When converting magnitudes to flux densities, the order of operations matters;
    the value of the unit needs to be established *before* the conversion. 
    For example, ``21 * u.ABmag.to(u.erg/u.s/u.cm**2/u.Hz)`` will give you 21 
    times :math:`f_\nu` for an AB mag of 1, whereas ``(21 * u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz)`` 
    will give you :math:`f_\nu` for an AB mag of 21.
 
Suppose we also knew the intrinsic color of the first star, then we can
calculate the reddening::

    >>> B_V0 = -0.2 * u.mag
    >>> EB_V = (B - V)[0] - B_V0
    >>> R_V = 3.1
    >>> A_V = R_V * EB_V
    >>> A_B = (R_V+1) * EB_V
    >>> EB_V, A_V, A_B  # doctest: +FLOAT_CMP
    (<Magnitude 0.4 mag>, <Quantity 1.24 mag>, <Quantity 1.64 mag>)

Here, one sees that the extinctions have been converted to quantities. This
happens generally for division and multiplication, since these processes
work only for dimensionless magnitudes (otherwise, the physical unit would have
to be raised to some power), and |quantity| objects, unlike logarithmic
quantities, allow units like ``mag / d``.

Note that one can take the automatic unit conversion quite far (perhaps too
far, but it is fun).  For instance, suppose we also knew the bolometric
correction and absolute bolometric magnitude, then we can calculate the
distance modulus::

    >>> BC_V = -0.3 * (u.m_bol - u.STmag)
    >>> M_bol = 5.46 * u.M_bol
    >>> DM = V[0] - A_V + BC_V - M_bol
    >>> BC_V, M_bol, DM  # doctest: +FLOAT_CMP
    (<Magnitude -0.3 mag(bol / ST)>,
     <Magnitude 5.46 mag(Bol)>,
     <Magnitude 10. mag(bol / Bol)>)

With a proper equivalency, we can also convert to distance without remembering
the 5-5log rule::

    >>> radius_and_inverse_area = [(u.pc, u.pc**-2,
    ...                            lambda x: 1./(4.*np.pi*x**2),
    ...                            lambda x: np.sqrt(1./(4.*np.pi*x)))]
    >>> DM.to(u.pc, equivalencies=radius_and_inverse_area)  # doctest: +FLOAT_CMP
    <Quantity 1000. pc>

Numpy functions
===============

For logarithmic quantities, most numpy functions and many array methods do not
make sense, hence they are disabled.  But one can use those one would expect to
work::

    >>> np.max(v_i)  # doctest: +FLOAT_CMP
    <Magnitude 4.00514998 mag(ct / s)>
    >>> np.std(v_i)  # doctest: +FLOAT_CMP
    <Magnitude 2.33971149 mag>

.. note:: This is implemented by having a list of supported ufuncs in
	  ``units/function/core.py`` and by explicitly disabling some
	  array methods in :class:`~astropy.units.function.FunctionQuantity`.
          If you believe a function or method is incorrectly treated,
	  please `let us know <http://www.astropy.org/contribute.html>`_.

Dimensionless logarithmic quantities
====================================

Dimensionless quantities are treated somewhat specially, in that, if needed,
logarithmic quantities will be converted to normal |quantity| objects with the
appropriate unit of ``mag``, ``dB``, or ``dex``.  With this, it is possible to
use composite units like ``mag/d`` or ``dB/m``, which cannot easily be
supported as logarithmic units.  For instance::

    >>> dBm = u.dB(u.mW)
    >>> signal_in, signal_out = 100. * dBm, 50 * dBm
    >>> cable_loss = (signal_in - signal_out) / (100. * u.m)
    >>> signal_in, signal_out, cable_loss  # doctest: +FLOAT_CMP
    (<Decibel 100. dB(mW)>, <Decibel 50. dB(mW)>, <Quantity 0.5 dB / m>)
    >>> better_cable_loss = 0.2 * u.dB / u.m
    >>> signal_in - better_cable_loss * 100. * u.m  # doctest: +FLOAT_CMP
    <Decibel 80. dB(mW)>


.. [M15] Mamajek et al., 2015, `arXiv:1510.06262
	  <https://ui.adsabs.harvard.edu/abs/2015arXiv151006262M>`_
.. [H95] E.g., Holtzman et al., 1995, `PASP 107, 1065
          <https://ui.adsabs.harvard.edu/abs/1995PASP..107.1065H>`_
.. [OG83] Oke, J.B., & Gunn, J. E., 1983, `ApJ 266, 713
	  <https://ui.adsabs.harvard.edu/abs/1983ApJ...266..713O>`_
