.. |quantity| replace:: :class:`~astropy.units.Quantity`

.. _unit_equivalencies:

Equivalencies
*************

The unit module has machinery for supporting equivalences between
different units in certain contexts, namely when equations can
uniquely relate a value in one unit to a different unit. A good
example is the equivalence between wavelength, frequency and energy
for specifying a wavelength of radiation. Normally these units are not
convertible, but when understood as representing light, they are
convertible in certain contexts.  Here we describe how to use the
equivalencies included in `astropy.units` and how to
define new equivalencies.

Equivalencies are used by passing a list of equivalency pairs to the
``equivalencies`` keyword argument of :meth:`Quantity.to
<astropy.units.quantity.Quantity.to>` or :meth:`Unit.to
<astropy.units.core.UnitBase.to>` methods.  Alternatively, if a larger
piece of code needs the same equivalencies, one can set them for a
:ref:`given context <equivalency-context>`.

Built-in equivalencies
======================

How to Convert Parallax to Distance
-----------------------------------

The length unit *parsec* is defined such that a star one parsec away
will exhibit a 1-arcsecond parallax. (Think of it as a contraction
between *parallax* and *arcsecond*.)

The :func:`~astropy.units.equivalencies.parallax` function handles
conversions between parallax angles and length.

In general, you should not be able to change units of length into
angles or vice versa, so :meth:`~astropy.units.core.UnitBase.to`
raises an exception::

  >>> from astropy import units as u
  >>> (0.8 * u.arcsec).to(u.parsec)  # doctest: +IGNORE_EXCEPTION_DETAIL
  Traceback (most recent call last):
    ...
  UnitConversionError: 'arcsec' (angle) and 'pc' (length) are not convertible

To trigger the conversion between parallax angle and distance, provide
:func:`~astropy.units.equivalencies.parallax` as the optional keyword
argument (``equivalencies=``) to the
:meth:`~astropy.units.core.UnitBase.to` method.

    >>> (0.8 * u.arcsec).to(u.parsec, equivalencies=u.parallax())
    <Quantity 1.25 pc>

Angles as Dimensionless Units
-----------------------------

Angles are treated as a physically distinct type, which usually helps
to avoid mistakes.  However, this is not very handy when working with
units related to rotational energy or the small angle approximation.
(Indeed, this double-sidedness underlies why radian went from
`supplementary to derived unit <https://www.bipm.org/en/CGPM/db/20/8/>`__.)
The function :func:`~astropy.units.equivalencies.dimensionless_angles`
provides the required equivalency list that helps convert between
angles and dimensionless units.  It is somewhat
different from all others in that it allows an arbitrary change in the
number of powers to which radian is raised (i.e., including zero and thus
dimensionless).  For instance, normally the following raise exceptions::

  >>> from astropy import units as u
  >>> u.degree.to('')  # doctest: +IGNORE_EXCEPTION_DETAIL
  Traceback (most recent call last):
    ...
  UnitConversionError: 'deg' (angle) and '' (dimensionless) are not convertible
  >>> (u.kg * u.m**2 * (u.cycle / u.s)**2).to(u.J)  # doctest: +IGNORE_EXCEPTION_DETAIL
  Traceback (most recent call last):
    ...
  UnitConversionError: 'cycle2 kg m2 / s2' and 'J' (energy) are not convertible

But when passing the proper conversion function,
:func:`~astropy.units.equivalencies.dimensionless_angles`, it works.

  >>> u.deg.to('', equivalencies=u.dimensionless_angles())  # doctest: +FLOAT_CMP
  0.017453292519943295
  >>> (0.5e38 * u.kg * u.m**2 * (u.cycle / u.s)**2).to(u.J,
  ...                            equivalencies=u.dimensionless_angles())  # doctest: +FLOAT_CMP
  <Quantity 1.9739208802178715e+39 J>
  >>> import numpy as np
  >>> np.exp((1j*0.125*u.cycle).to('', equivalencies=u.dimensionless_angles()))  # doctest: +SKIP
  <Quantity  0.70710678+0.70710678j>

The example with complex numbers is also one may well be doing a fair
number of similar calculations.  For such situations, there is the
option to :ref:`set default equivalencies <equivalency-context>`.

In some situations, this equivalency may behave differently than
anticipated.  For instance, it might at first seem reasonable to use it
for converting from an angular velocity :math:`\omega` in radians per
second to the corresponding frequency :math:`f` in hertz (i.e., to
implement :math:`f=\omega/2\pi`). However, attempting this yields:

  >>> (1*u.rad/u.s).to(u.Hz, equivalencies=u.dimensionless_angles())  # doctest: +FLOAT_CMP
  <Quantity 1. Hz>
  >>> (1*u.cycle/u.s).to(u.Hz, equivalencies=u.dimensionless_angles())  # doctest: +FLOAT_CMP
  <Quantity 6.283185307179586 Hz>

Here, we might have expected ~0.159 Hz in the first example and 1 Hz in
the second. However, :func:`~astropy.units.equivalencies.dimensionless_angles`
converts to radians per second and then drops radians as a unit. The
implicit mistake made in these examples is that the unit Hz is taken to be
equivalent to cycles per second, which it is not (it is just "per second").
This realization also leads to the solution: to use an explicit equivalency
between cycles per second and hertz:

  >>> (1*u.rad/u.s).to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])  # doctest: +FLOAT_CMP
  <Quantity 0.15915494309189535 Hz>
  >>> (1*u.cy/u.s).to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])  # doctest: +FLOAT_CMP
  <Quantity 1. Hz>

Spectral Units
--------------

:func:`~astropy.units.equivalencies.spectral` is a function that returns
an equivalency list to handle conversions between wavelength,
frequency, energy, and wave number.

As mentioned above with parallax units, we simply pass a list of
equivalencies (in this case, the result of
:func:`~astropy.units.equivalencies.spectral`) as the third argument to the
:meth:`~astropy.units.core.UnitBase.to` method and wavelength, frequency and
energy can be converted.

  >>> ([1000, 2000] * u.nm).to(u.Hz, equivalencies=u.spectral())  # doctest: +FLOAT_CMP
  <Quantity [2.99792458e+14, 1.49896229e+14] Hz>
  >>> ([1000, 2000] * u.nm).to(u.eV, equivalencies=u.spectral())  # doctest: +FLOAT_CMP
  <Quantity [1.23984193, 0.61992096] eV>

These equivalencies even work with non-base units::

  >>> # Inches to calories
  >>> from astropy.units import imperial
  >>> imperial.inch.to(imperial.Cal, equivalencies=u.spectral())  # doctest: +FLOAT_CMP
  1.869180759162485e-27

.. _astropy-units-dopper-equivalencies:

Spectral (Doppler) equivalencies
--------------------------------

Spectral equivalencies allow you to convert between wavelength,
frequency, energy, and wave number but not to velocity, which is
frequently the quantity of interest.

It is fairly straightforward to define the equivalency, but note that there are
different `conventions <https://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`__.
In these conventions :math:`f_0` is the rest frequency, :math:`f` is the observed frequency,
:math:`V` is the velocity, and :math:`c` is the speed of light:

    * Radio         :math:`V = c \frac{f_0 - f}{f_0}  ;  f(V) = f_0 ( 1 - V/c )`
    * Optical       :math:`V = c \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + V/c )^{-1}`
    * Relativistic  :math:`V = c \frac{f_0^2 - f^2}{f_0^2 + f^2} ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

These three conventions are implemented in
:mod:`astropy.units.equivalencies` as
:func:`~astropy.units.equivalencies.doppler_optical`,
:func:`~astropy.units.equivalencies.doppler_radio`, and
:func:`~astropy.units.equivalencies.doppler_relativistic`.  Example use::

    >>> restfreq = 115.27120 * u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> freq_to_vel = u.doppler_radio(restfreq)
    >>> (116e9 * u.Hz).to(u.km / u.s, equivalencies=freq_to_vel)  # doctest: +FLOAT_CMP
    <Quantity -1895.4321928669085 km / s>

Spectral Flux / Luminosity Density Units
----------------------------------------

There is also support for spectral flux and luminosity density units
and their equivalent surface brightness units. Their use
is more complex, since it is necessary to also supply the location in the
spectrum for which the conversions will be done, and the units of those spectral
locations.  The function that handles these unit conversions is
:func:`~astropy.units.equivalencies.spectral_density`. This function takes as
its arguments the |quantity| for the spectral location. For example::

    >>> (1.5 * u.Jy).to(u.photon / u.cm**2 / u.s / u.Hz,
    ...                 equivalencies=u.spectral_density(3500 * u.AA)) # doctest: +FLOAT_CMP
    <Quantity 2.6429114293019694e-12 ph / (cm2 Hz s)>
    >>> (1.5 * u.Jy).to(u.photon / u.cm**2 / u.s / u.micron,
    ...                 equivalencies=u.spectral_density(3500 * u.AA))  # doctest: +FLOAT_CMP
    <Quantity 6467.9584789120845 ph / (cm2 micron s)>
    >>> a = 1. * u.photon / u.s / u.angstrom
    >>> a.to(u.erg / u.s / u.Hz,
    ...      equivalencies=u.spectral_density(5500 * u.AA)) # doctest: +FLOAT_CMP
    <Quantity 3.6443382634999996e-23 erg / (Hz s)>

Brightness Temperature / Surface Brightness Equivalency
-------------------------------------------------------

There is an equivalency between surface brightness (flux density per area) and
brightness temperature.  This equivalency is often referred to as "Antenna
Gain" since, at a given frequency, telescope brightness sensitivity is
unrelated to aperture size, but flux density sensitivity is, so this
equivalency is only dependent on the aperture size.  See `Tools of Radio
Astronomy
<https://books.google.com/books?id=9KHw6R8rQEMC&pg=PA179&source=gbs_toc_r&cad=4#v=onepage&q&f=false>`__
for details.

.. note:: The brightness temperature mentioned here is the Rayleigh-Jeans
          equivalent temperature, which results in a linear relation between
          flux and temperature. This is the convention that is most often used
          in relation to observations, but if you are interested in computing
          the *exact* temperature of a planck function that would produce a
          given flux, you should not use this equivalency.

The `~astropy.units.equivalencies.brightness_temperature` equivalency requires
the beam area and frequency as arguments.  Recalling that the area of a 2D
gaussian is :math:`2 \pi \sigma^2` (see `wikipedia
<https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function>`_),
here is an example::

    >>> import numpy as np
    >>> beam_sigma = 50*u.arcsec
    >>> omega_B = 2 * np.pi * beam_sigma**2
    >>> freq = 5 * u.GHz
    >>> (1*u.Jy/omega_B).to(u.K, equivalencies=u.brightness_temperature(freq))  # doctest: +FLOAT_CMP
    <Quantity 3.526295144567176 K>

If you have beam full-width half-maxima (FWHM), which are often quoted and are
the values stored in the FITS header keywords BMAJ and BMIN, a more appropriate
example converts the FWHM to sigma::

    >>> import numpy as np
    >>> beam_fwhm = 50*u.arcsec
    >>> fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
    >>> beam_sigma = beam_fwhm * fwhm_to_sigma
    >>> omega_B = 2 * np.pi * beam_sigma**2
    >>> freq = 5 * u.GHz
    >>> (1*u.Jy/omega_B).to(u.K, equivalencies=u.brightness_temperature(freq))  # doctest: +FLOAT_CMP
    <Quantity 19.553932298231704 K>

You can also convert between ``Jy/beam`` and ``K`` by specifying the beam area::

    >>> import numpy as np
    >>> beam_fwhm = 50*u.arcsec
    >>> fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
    >>> beam_sigma = beam_fwhm * fwhm_to_sigma
    >>> omega_B = 2 * np.pi * beam_sigma**2
    >>> freq = 5 * u.GHz
    >>> (1*u.Jy/u.beam).to(u.K, u.brightness_temperature(freq, beam_area=omega_B))  # doctest: +FLOAT_CMP
    <Quantity 19.553932298231704 K>

Beam Equivalency
----------------

Radio data, especially from interferometers, is often produced in units of ``Jy/beam``.
Converting this number to a beam-independent value, e.g., ``Jy/sr``, can be done
with the `~astropy.units.equivalencies.beam_angular_area` equivalency::

    >>> import numpy as np
    >>> beam_fwhm = 50*u.arcsec
    >>> fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
    >>> beam_sigma = beam_fwhm * fwhm_to_sigma
    >>> omega_B = 2 * np.pi * beam_sigma**2
    >>> (1*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(omega_B))  # doctest: +FLOAT_CMP
    <Quantity 15.019166691021288 MJy / sr>


Note that the `radio_beam <https://github.com/radio-astro-tools/radio-beam>`_ package
deals with beam input/output and various operations more directly.

Temperature Energy Equivalency
------------------------------

This equivalency allows conversion between temperature and its equivalent
in energy (i.e., the temperature multiplied by the Boltzmann constant),
usually expressed in electronvolts. This is used frequently for
observations at high-energy, be it for solar or X-ray astronomy. Example::

    >>> import astropy.units as u
    >>> t_k = 1e6 * u.K
    >>> t_k.to(u.eV, equivalencies=u.temperature_energy())  # doctest: +FLOAT_CMP
    <Quantity 86.17332384960955 eV>

.. _tcmb-equivalency:

Thermodynamic Temperature Equivalency
-------------------------------------

This :func:`~astropy.units.equivalencies.thermodynamic_temperature`
equivalency allows conversion between Jy/beam and "thermodynamic
temperature", :math:`T_{CMB}`, in Kelvins. Example::

    >>> import astropy.units as u
    >>> nu = 143 * u.GHz
    >>> t_k = 0.002632051878 * u.K
    >>> t_k.to(u.MJy / u.sr, equivalencies=u.thermodynamic_temperature(nu))  # doctest: +FLOAT_CMP
    <Quantity 1. MJy / sr>

By default, this will use the :math:`T_{CMB}` value for the 'default cosmology'
in astropy, but it is possible to specify a custom :math:`T_{CMB}` value for a
specific cosmology as the second argument to the equivalency::

    >>> from astropy.cosmology import WMAP9
    >>> t_k.to(u.MJy / u.sr, equivalencies=u.thermodynamic_temperature(nu, T_cmb=WMAP9.Tcmb0))  # doctest: +FLOAT_CMP
    <Quantity 0.99982392 MJy / sr>

Molar Mass AMU Equivalency
--------------------------

This equivalency allows conversion
between the atomic mass unit and the equivalent g/mol.
For reference to why this was added,
refer to `astropy GitHub issue 6040 <https://github.com/astropy/astropy/issues/6040>`_
The following is an example of it's usage:

    >>> import astropy.units as u
    >>> import astropy.constants as const
    >>> x = 1 * (u.g / u.mol)
    >>> y = 1 * u.u
    >>> x.to(u.u, equivalencies=u.molar_mass_amu()) # doctest: +FLOAT_CMP
    <Quantity 1.0 u>
    >>> y.to(u.g/u.mol, equivalencies=u.molar_mass_amu()) # doctest: +FLOAT_CMP
    <Quantity 1.0 g / mol>

Pixel and plate scale Equivalencies
-----------------------------------

These equivalencies are for converting between angular scales and either linear
scales in the focal plane or distances in units of the number of pixels.  For
example, suppose you are working with cutouts from the Sloan Digital Sky Survey,
which defaults to a pixel scale of 0.4 arcseconds per pixel, and want to know
the true size of something that you measure to be 240 pixels across in the
cutout image::

    >>> import astropy.units as u
    >>> sdss_pixelscale = u.pixel_scale(0.4*u.arcsec/u.pixel)
    >>> (240*u.pixel).to(u.arcmin, sdss_pixelscale)  # doctest: +FLOAT_CMP
    <Quantity 1.6 arcmin>

Or maybe you are designing an instrument for a telescope that someone told you
has a (inverse) plate  scale of 7.8 meters per radian (for your desired focus),
and you want to know how big your pixels need to be to cover half an arcsecond::

    >>> import astropy.units as u
    >>> tel_platescale = u.plate_scale(7.8*u.m/u.radian)
    >>> (0.5*u.arcsec).to(u.micron, tel_platescale)  # doctest: +FLOAT_CMP
    <Quantity 18.9077335632719 micron>

Photometric Zero Point Equivalency
----------------------------------

This equivalency provides an easy way to move between photometric systems (i.e.,
those defined relative to a particular zero-point flux) and absolute fluxes.
This is most useful in conjuction with support for :ref:`logarithmic_units`.
For example, suppose you are observing a target with a filter with a reported
standard zero point of 3631.1 Jy::

    >>> target_flux = 1.2 * u.nanomaggy
    >>> zero_point_star_equiv = u.zero_point_flux(3631.1 * u.Jy)
    >>> u.Magnitude(target_flux.to(u.AB, zero_point_star_equiv))  # doctest: +FLOAT_CMP
    <Magnitude 22.30195136 mag(AB)>

.. _H0-equivalency:

Reduced Hubble constant/"little-h" Equivalency
----------------------------------------------

The dimensionless version of the Hubble constant - often known as "little h" -
is a frequently-used quantity in extragalactic astrophysics.  It is also widely
known as the bane of beginners' existence in such fields (See e.g., the title of
`this paper <https://doi.org/10.1017/pasa.2013.31>`__, which also provides
valuable advice on the use of little h).  Astropy provides an equivalency that
helps keep this straight in at least some of these cases, by providing a way to
convert to/from physical to "little h" units.  Two example conversions:

    >>> import astropy.units as u
    >>> H0_70 = 70 * u.km/u.s / u.Mpc
    >>> distance = 70 * (u.Mpc/u.littleh)
    >>> distance.to(u.Mpc, u.with_H0(H0_70))  # doctest: +FLOAT_CMP
    <Quantity 100.0 Mpc>
    >>> luminosity = 0.49 * u.Lsun * u.littleh**-2
    >>> luminosity.to(u.Lsun, u.with_H0(H0_70))  # doctest: +FLOAT_CMP
    <Quantity 1.0 solLum>

Note the unit name ``littleh`` - while this unit is usually expressed in the
literature as just ``h``, here it is ``littleh`` to not cause confusion with
"hours".

If no argument is given (or the argument is `None`), this equivalency assumes
the ``H0`` from the current default cosmology:

    >>> distance = 100 * (u.Mpc/u.littleh)
    >>> distance.to(u.Mpc, u.with_H0())  # doctest: +FLOAT_CMP
    <Quantity 147.62326543 Mpc>

This equivalency also allows a common magnitude formulation of little h
scaling:

    >>> mag_quantity = 12 * (u.mag - u.MagUnit(u.littleh**2))
    >>> mag_quantity  # doctest: +FLOAT_CMP
    <Magnitude 12. mag(1 / littleh2)>
    >>> mag_quantity.to(u.mag, u.with_H0(H0_70))  # doctest: +FLOAT_CMP
    <Quantity 11.2254902 mag>

Temperature Equivalency
-----------------------

The :func:`~astropy.units.temperature` equivalency allows conversion
between the Celsius, Fahrenheit, and Kelvin. For example::

    >>> import astropy.units as u
    >>> temp_C = 0 * u.Celsius
    >>> temp_Kelvin = temp_C.to(u.K, equivalencies=u.temperature())
    >>> temp_Kelvin  # doctest: +FLOAT_CMP
    <Quantity 273.15 K>
    >>> temp_F = temp_C.to(u.imperial.deg_F, equivalencies=u.temperature())
    >>> temp_F  # doctest: +FLOAT_CMP
    <Quantity 32. deg_F>

.. note:: You can also use ``u.deg_C`` instead of ``u.Celsius``.

Mass-Energy Equivalency
-----------------------

In a special relativity context, mass and energy can be equivalent units. For
example::

    >>> import astropy.units as u
    >>> (1 * u.g).to(u.eV, u.mass_energy())  # doctest: +FLOAT_CMP
    <Quantity 5.60958865e+32 eV>


Writing new equivalencies
=========================

An equivalence list is just a list of tuples, where each tuple has 4
elements::

  (from_unit, to_unit, forward, backward)

``from_unit`` and ``to_unit`` are the equivalent units.  ``forward`` and
``backward`` are functions that convert values between those units. ``forward``
and ``backward`` are optional, and if omitted such an equivalency simply
declares that the two units should be taken as equivalent.

For example, until 1964 the metric liter was defined as the volume of
1kg of water at 4°C at 760mm mercury pressure.  Volumes and masses are
not normally directly convertible, but if we hold the constants in the
1964 definition of the liter as true, we could build an equivalency
for them::

  >>> liters_water = [
  ...    (u.l, u.g, lambda x: 1000.0 * x, lambda x: x / 1000.0)
  ... ]
  >>> u.l.to(u.kg, 1, equivalencies=liters_water)
  1.0

Note that the equivalency can be used with any other compatible units::

  >>> from astropy.units import imperial
  >>> imperial.gallon.to(imperial.pound, 1, equivalencies=liters_water)  # doctest: +FLOAT_CMP
  8.345404463333525

And it also works in the other direction::

  >>> imperial.lb.to(imperial.pint, 1, equivalencies=liters_water)  # doctest: +FLOAT_CMP
  0.9586114172355459

A slightly more complicated example: Spectral Doppler Equivalencies
-------------------------------------------------------------------

We show how to define an equivalency using the radio convention for CO 1-0.
This function is already defined in
:func:`~astropy.units.equivalencies.doppler_radio`,
but this example is illustrative::

    >>> from astropy.constants import si
    >>> restfreq = 115.27120  # rest frequency of 12 CO 1-0 in GHz
    >>> freq_to_vel = [(u.GHz, u.km/u.s,
    ... lambda x: (restfreq-x) / restfreq * si.c.to_value('km/s'),
    ... lambda x: (1-x/si.c.to_value('km/s')) * restfreq )]
    >>> u.Hz.to(u.km / u.s, 116e9, equivalencies=freq_to_vel)  # doctest: +FLOAT_CMP
    -1895.4321928669262
    >>> (116e9 * u.Hz).to(u.km / u.s, equivalencies=freq_to_vel)  # doctest: +FLOAT_CMP
    <Quantity -1895.4321928669262 km / s>

Note that once this is defined for GHz and km/s, it will work for all other
units of frequency and velocity.  ``x`` is converted from the input frequency
unit (e.g., Hz) to GHz before being passed to ``lambda x:``.  Similarly, the
return value is assumed to be in units of ``km/s``, which is why the ``.value``
of ``c`` is used instead of the constant.

Displaying available equivalencies
==================================

The :meth:`~astropy.units.core.UnitBase.find_equivalent_units` method also
understands equivalencies.  For example, without passing equivalencies,
there are three compatible units for ``Hz`` in the standard set::

  >>> u.Hz.find_equivalent_units()
    Primary name | Unit definition | Aliases
  [
    Bq           | 1 / s           | becquerel    ,
    Ci           | 3.7e+10 / s    | curie        ,
    Hz           | 1 / s           | Hertz, hertz ,
  ]

However, when passing the spectral equivalency, you can see there are
all kinds of things that ``Hz`` can be converted to::

  >>> u.Hz.find_equivalent_units(equivalencies=u.spectral())
    Primary name | Unit definition        | Aliases
  [
    AU           | 1.49598e+11 m          | au, astronomical_unit ,
    Angstrom     | 1e-10 m                | AA, angstrom          ,
    Bq           | 1 / s                  | becquerel             ,
    Ci           | 3.7e+10 / s            | curie                 ,
    Hz           | 1 / s                  | Hertz, hertz          ,
    J            | kg m2 / s2             | Joule, joule          ,
    Ry           | 2.17987e-18 kg m2 / s2 | rydberg               ,
    cm           | 0.01 m                 | centimeter            ,
    eV           | 1.60218e-19 kg m2 / s2 | electronvolt          ,
    earthRad     | 6.3781e+06 m           | R_earth, Rearth       ,
    erg          | 1e-07 kg m2 / s2       |                       ,
    jupiterRad   | 7.1492e+07 m           | R_jup, Rjup, R_jupiter, Rjupiter ,
    k            | 100 / m                | Kayser, kayser        ,
    lyr          | 9.46073e+15 m          | lightyear             ,
    m            | irreducible            | meter                 ,
    micron       | 1e-06 m                |                       ,
    pc           | 3.08568e+16 m          | parsec                ,
    solRad       | 6.957e+08 m            | R_sun, Rsun           ,
  ]

.. _equivalency-context:

Using equivalencies in larger pieces of code
============================================
Sometimes one has an involved calculation where one is regularly
switching back between equivalent units. For these cases, one can set
equivalencies that will by default be used, in a way similar to which
one can :ref:`enable other units <enabling-other-units>`.

For instance, to enable radian to be treated as a dimensionless unit,
simply do:

.. doctest-skip::

  >>> import astropy.units as u
  >>> u.set_enabled_equivalencies(u.dimensionless_angles())
  <astropy.units.core._UnitContext object at ...>
  >>> u.deg.to('')  # doctest: +FLOAT_CMP
  0.017453292519943295

Here, any list of equivalencies could be used, or one could add, e.g.,
:func:`~astropy.units.equivalencies.spectral` and
:func:`~astropy.units.equivalencies.spectral_density` (since these return
lists, they should indeed be combined by adding them together).

The disadvantage of the above approach is that you may forget to turn
the default off (done by giving an empty argument). To automate this,
a context manager is provided:

.. doctest-skip::

  >>> import astropy.units as u
  >>> with u.set_enabled_equivalencies(u.dimensionless_angles()):
  ...    phase = 0.5 * u.cycle
  ...    c = np.exp(1j*phase)
  >>> c  # doctest: +FLOAT_CMP
  <Quantity (-1+1.2246063538223773e-16j) >
