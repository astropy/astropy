# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Models that have physical origins.
"""
# pylint: disable=invalid-name, no-member

import warnings

import numpy as np

from astropy import constants as const
from astropy import units as u
from astropy import cosmology
from astropy.utils.exceptions import AstropyUserWarning
from .core import Fittable1DModel
from .parameters import Parameter, InputParameterError

__all__ = ["BlackBody", "Drude1D", "Plummer1D", "NFW"]


class BlackBody(Fittable1DModel):
    """
    Blackbody model using the Planck function.

    Parameters
    ----------
    temperature : :class:`~astropy.units.Quantity`
        Blackbody temperature.

    scale : float or :class:`~astropy.units.Quantity`
        Scale factor

    Notes
    -----

    Model formula:

        .. math:: B_{\\nu}(T) = A \\frac{2 h \\nu^{3} / c^{2}}{exp(h \\nu / k T) - 1}

    Examples
    --------
    >>> from astropy.modeling import models
    >>> from astropy import units as u
    >>> bb = models.BlackBody(temperature=5000*u.K)
    >>> bb(6000 * u.AA)  # doctest: +FLOAT_CMP
    <Quantity 1.53254685e-05 erg / (cm2 Hz s sr)>

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import BlackBody
        from astropy import units as u
        from astropy.visualization import quantity_support

        bb = BlackBody(temperature=5778*u.K)
        wav = np.arange(1000, 110000) * u.AA
        flux = bb(wav)

        with quantity_support():
            plt.figure()
            plt.semilogx(wav, flux)
            plt.axvline(bb.nu_max.to(u.AA, equivalencies=u.spectral()).value, ls='--')
            plt.show()
    """

    # We parametrize this model with a temperature and a scale.
    temperature = Parameter(default=5000.0, min=0, unit=u.K)
    scale = Parameter(default=1.0, min=0)

    # We allow values without units to be passed when evaluating the model, and
    # in this case the input x values are assumed to be frequencies in Hz.
    _input_units_allow_dimensionless = True

    # We enable the spectral equivalency by default for the spectral axis
    input_units_equivalencies = {'x': u.spectral()}

    def evaluate(self, x, temperature, scale):
        """Evaluate the model.

        Parameters
        ----------
        x : float, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Frequency at which to compute the blackbody. If no units are given,
            this defaults to Hz.

        temperature : float, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Temperature of the blackbody. If no units are given, this defaults
            to Kelvin.

        scale : float, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Desired scale for the blackbody.

        Returns
        -------
        y : number or ndarray
            Blackbody spectrum. The units are determined from the units of
            ``scale``.

        .. note::

            Use `numpy.errstate` to suppress Numpy warnings, if desired.

        .. warning::

            Output values might contain ``nan`` and ``inf``.

        Raises
        ------
        ValueError
            Invalid temperature.

        ZeroDivisionError
            Wavelength is zero (when converting to frequency).
        """
        if not isinstance(temperature, u.Quantity):
            in_temp = u.Quantity(temperature, u.K)
        else:
            in_temp = temperature

        # Convert to units for calculations, also force double precision
        with u.add_enabled_equivalencies(u.spectral() + u.temperature()):
            freq = u.Quantity(x, u.Hz, dtype=np.float64)
            temp = u.Quantity(in_temp, u.K)

        # check the units of scale and setup the output units
        bb_unit = u.erg / (u.cm ** 2 * u.s * u.Hz * u.sr)  # default unit
        # use the scale that was used at initialization for determining the units to return
        # to support returning the right units when fitting where units are stripped
        if hasattr(self.scale, "unit") and self.scale.unit is not None:
            # check that the units on scale are covertable to surface brightness units
            if not self.scale.unit.is_equivalent(bb_unit, u.spectral_density(x)):
                raise ValueError(
                    f"scale units not surface brightness: {self.scale.unit}"
                )
            # use the scale passed to get the value for scaling
            if hasattr(scale, "unit"):
                mult_scale = scale.value
            else:
                mult_scale = scale
            bb_unit = self.scale.unit
        else:
            mult_scale = scale

        # Check if input values are physically possible
        if np.any(temp < 0):
            raise ValueError(f"Temperature should be positive: {temp}")
        if not np.all(np.isfinite(freq)) or np.any(freq <= 0):
            warnings.warn(
                "Input contains invalid wavelength/frequency value(s)",
                AstropyUserWarning,
            )

        log_boltz = const.h * freq / (const.k_B * temp)
        boltzm1 = np.expm1(log_boltz)

        # Calculate blackbody flux
        bb_nu = 2.0 * const.h * freq ** 3 / (const.c ** 2 * boltzm1) / u.sr

        y = mult_scale * bb_nu.to(bb_unit, u.spectral_density(freq))

        # If the temperature parameter has no unit, we should return a unitless
        # value. This occurs for instance during fitting, since we drop the
        # units temporarily.
        if hasattr(temperature, "unit"):
            return y
        return y.value

    @property
    def input_units(self):
        # The input units are those of the 'x' value, which should always be
        # Hz. Because we do this, and because input_units_allow_dimensionless
        # is set to True, dimensionless values are assumed to be in Hz.
        return {self.inputs[0]: u.Hz}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {"temperature": u.K}

    @property
    def bolometric_flux(self):
        """Bolometric flux."""
        # bolometric flux in the native units of the planck function
        native_bolflux = (
            self.scale.value * const.sigma_sb * self.temperature ** 4 / np.pi
        )
        # return in more "astro" units
        return native_bolflux.to(u.erg / (u.cm ** 2 * u.s))

    @property
    def lambda_max(self):
        """Peak wavelength when the curve is expressed as power density."""
        return const.b_wien / self.temperature

    @property
    def nu_max(self):
        """Peak frequency when the curve is expressed as power density."""
        return 2.8214391 * const.k_B * self.temperature / const.h


class Drude1D(Fittable1DModel):
    """
    Drude model based one the behavior of electons in materials (esp. metals).

    Parameters
    ----------
    amplitude : float
        Peak value
    x_0 : float
        Position of the peak
    fwhm : float
        Full width at half maximum

    Model formula:

        .. math:: f(x) = A \\frac{(fwhm/x_0)^2}{((x/x_0 - x_0/x)^2 + (fwhm/x_0)^2}

    Examples
    --------

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Drude1D

        fig, ax = plt.subplots()

        # generate the curves and plot them
        x = np.arange(7.5 , 12.5 , 0.1)

        dmodel = Drude1D(amplitude=1.0, fwhm=1.0, x_0=10.0)
        ax.plot(x, dmodel(x))

        ax.set_xlabel('x')
        ax.set_ylabel('F(x)')

        plt.show()
    """

    amplitude = Parameter(default=1.0)
    x_0 = Parameter(default=1.0)
    fwhm = Parameter(default=1.0)

    @staticmethod
    def evaluate(x, amplitude, x_0, fwhm):
        """
        One dimensional Drude model function
        """
        return (
            amplitude
            * ((fwhm / x_0) ** 2)
            / ((x / x_0 - x_0 / x) ** 2 + (fwhm / x_0) ** 2)
        )

    @staticmethod
    def fit_deriv(x, amplitude, x_0, fwhm):
        """
        Drude1D model function derivatives.
        """
        d_amplitude = (fwhm / x_0) ** 2 / ((x / x_0 - x_0 / x) ** 2 + (fwhm / x_0) ** 2)
        d_x_0 = (
            -2
            * amplitude
            * d_amplitude
            * (
                (1 / x_0)
                + d_amplitude
                * (x_0 ** 2 / fwhm ** 2)
                * (
                    (-x / x_0 - 1 / x) * (x / x_0 - x_0 / x)
                    - (2 * fwhm ** 2 / x_0 ** 3)
                )
            )
        )
        d_fwhm = (2 * amplitude * d_amplitude / fwhm) * (1 - d_amplitude)
        return [d_amplitude, d_x_0, d_fwhm]

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "fwhm": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }

    @property
    def return_units(self):
        if self.amplitude.unit is None:
            return None
        return {self.outputs[0]: self.amplitude.unit}

    @x_0.validator
    def x_0(self, val):
        """ Ensure `x_0` is not 0."""
        if val == 0:
            raise InputParameterError("0 is not an allowed value for x_0")

    def bounding_box(self, factor=50):
        """Tuple defining the default ``bounding_box`` limits,
        ``(x_low, x_high)``.

        Parameters
        ----------
        factor : float
            The multiple of FWHM used to define the limits.
        """
        x0 = self.x_0
        dx = factor * self.fwhm

        return (x0 - dx, x0 + dx)


class Plummer1D(Fittable1DModel):
    r"""One dimensional Plummer density profile model.

    Parameters
    ----------
    mass : float
        Total mass of cluster.
    r_plum : float
        Scale parameter which sets the size of the cluster core.

    Notes
    -----
    Model formula:

    .. math::

        \rho(r)=\frac{3M}{4\pi a^3}(1+\frac{r^2}{a^2})^{-5/2}

    References
    ----------
    .. [1] https://ui.adsabs.harvard.edu/abs/1911MNRAS..71..460P
    """

    mass = Parameter(default=1.0)
    r_plum = Parameter(default=1.0)

    @staticmethod
    def evaluate(x, mass, r_plum):
        """
        Evaluate plummer density profile model.
        """
        return (3*mass)/(4 * np.pi * r_plum**3) * (1+(x/r_plum)**2)**(-5/2)

    @staticmethod
    def fit_deriv(x, mass, r_plum):
        """
        Plummer1D model derivatives.
        """
        d_mass = 3 / ((4*np.pi*r_plum**3) * (((x/r_plum)**2 + 1)**(5/2)))
        d_r_plum = (6*mass*x**2-9*mass*r_plum**2) / ((4*np.pi * r_plum**6) *
                                                     (1+(x/r_plum)**2)**(7/2))
        return [d_mass, d_r_plum]

    @property
    def input_units(self):
        if self.mass.unit is None and self.r_plum.unit is None:
            return None
        else:
            return {self.inputs[0]: self.r_plum.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'mass': outputs_unit[self.outputs[0]] * inputs_unit[self.inputs[0]] ** 3,
                'r_plum': inputs_unit[self.inputs[0]]}


class NFW(Fittable1DModel):
    r"""
    Navarro–Frenk–White (NFW) profile - model for radial distribution of dark matter.

    Parameters
    ----------
    mass : float or :class:`~astropy.units.Quantity`
        Mass of NFW peak within specified overdensity radius.
    concentration : float
        Concentration of the NFW profile.
    redshift : float
        Redshift of the NFW profile.
    massfactor : tuple or str
        Mass overdensity factor and type for provided profiles:
            Tuple version:
                ("virial",) : virial radius

                ("critical", N)  : radius where density is N times that of the critical density

                ("mean", N)  : radius where density is N times that of the mean density

            String version:
                "virial" : virial radius

                "Nc"  : radius where density is N times that of the critical density (e.g. "200c")

                "Nm"  : radius where density is N times that of the mean density (e.g. "500m")
    cosmo : :class:`~astropy.cosmology.Cosmology`
        Background cosmology for density calculation. If None, the default cosmology will be used.

    Notes
    -----

    Model formula:

    .. math:: \rho(r)=\frac{\delta_c\rho_{c}}{r/r_s(1+r/r_s)^2}

    References
    ----------
    .. [1] https://arxiv.org/pdf/astro-ph/9508025
    .. [2] https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile
    .. [3] https://en.wikipedia.org/wiki/Virial_mass
    """

    # Model Parameters

    # NFW Profile mass
    mass = Parameter(default=1.0, min=1.0, unit=u.M_sun)

    # NFW profile concentration
    concentration = Parameter(default=1.0, min=1.0)

    # NFW Profile redshift
    redshift = Parameter(default=0.0, min=0.0)

    # We allow values without units to be passed when evaluating the model, and
    # in this case the input r values are assumed to be lengths / positions in kpc.
    _input_units_allow_dimensionless = True

    def __init__(self, mass=u.Quantity(mass.default, mass.unit),
                 concentration=concentration.default, redshift=redshift.default,
                 massfactor=("critical", 200), cosmo=None,  **kwargs):
        # Set default cosmology
        if cosmo is None:
            cosmo = cosmology.default_cosmology.get()

        # Set mass overdensity type and factor
        self._density_delta(massfactor, cosmo, redshift)

        # Establish mass units for density calculation (default solar masses)
        if not isinstance(mass, u.Quantity):
            in_mass = u.Quantity(mass, u.M_sun)
        else:
            in_mass = mass

        # Obtain scale radius
        self._radius_s(mass, concentration)

        # Obtain scale density
        self._density_s(mass, concentration)

        super().__init__(mass=in_mass, concentration=concentration, redshift=redshift, **kwargs)

    def evaluate(self, r, mass, concentration, redshift):
        """
        One dimensional NFW profile function

        Parameters
        ----------
        r : float or :class:`~astropy.units.Quantity`
            Radial position of density to be calculated for the NFW profile.
        mass : float or :class:`~astropy.units.Quantity`
            Mass of NFW peak within specified overdensity radius.
        concentration : float
            Concentration of the NFW profile.
        redshift : float
            Redshift of the NFW profile.

        Returns
        -------
        density : float or :class:`~astropy.units.Quantity`
            NFW profile mass density at location ``r``. The density units are:
            [``mass`` / ``r`` ^3]

        Notes
        -----
        .. warning::

            Output values might contain ``nan`` and ``inf``.
        """
        # Create radial version of input with dimension
        if hasattr(r, "unit"):
            in_r = r
        else:
            in_r = u.Quantity(r, u.kpc)

        # Define reduced radius (r / r_{\\rm s})
        #   also update scale radius
        radius_reduced = in_r / self._radius_s(mass, concentration).to(in_r.unit)

        # Density distribution
        # \rho (r)=\frac{\rho_0}{\frac{r}{R_s}\left(1~+~\frac{r}{R_s}\right)^2}
        #   also update scale density
        density = self._density_s(mass, concentration) / (radius_reduced *
                                                          (u.Quantity(1.0) + radius_reduced) ** 2)

        if hasattr(mass, "unit"):
            return density
        else:
            return density.value

    def _density_delta(self, massfactor, cosmo, redshift):
        """
        Calculate density delta.
        """
        # Set mass overdensity type and factor
        if isinstance(massfactor, tuple):
            # Tuple options
            #   ("virial")       : virial radius
            #   ("critical", N)  : radius where density is N that of the critical density
            #   ("mean", N)      : radius where density is N that of the mean density
            if massfactor[0].lower() == "virial":
                # Virial Mass
                delta = None
                masstype = massfactor[0].lower()
            elif massfactor[0].lower() == "critical":
                # Critical or Mean Overdensity Mass
                delta = float(massfactor[1])
                masstype = 'c'
            elif massfactor[0].lower() == "mean":
                # Critical or Mean Overdensity Mass
                delta = float(massfactor[1])
                masstype = 'm'
            else:
                raise ValueError("Massfactor '" + str(massfactor[0]) + "' not one of 'critical', "
                                                                       "'mean', or 'virial'")
        else:
            try:
                # String options
                #   virial : virial radius
                #   Nc  : radius where density is N that of the critical density
                #   Nm  : radius where density is N that of the mean density
                if massfactor.lower() == "virial":
                    # Virial Mass
                    delta = None
                    masstype = massfactor.lower()
                elif massfactor[-1].lower() == 'c' or massfactor[-1].lower() == 'm':
                    # Critical or Mean Overdensity Mass
                    delta = float(massfactor[0:-1])
                    masstype = massfactor[-1].lower()
                else:
                    raise ValueError("Massfactor " + str(massfactor) + " string not of the form "
                                                                       "'#m', '#c', or 'virial'")
            except (AttributeError, TypeError):
                raise TypeError("Massfactor " + str(
                    massfactor) + " not a tuple or string")

        # Set density from masstype specification
        if masstype == "virial":
            Om_c = cosmo.Om(redshift) - 1.0
            d_c = 18.0 * np.pi ** 2 + 82.0 * Om_c - 39.0 * Om_c ** 2
            self.density_delta = d_c * cosmo.critical_density(redshift)
        elif masstype == 'c':
            self.density_delta = delta * cosmo.critical_density(redshift)
        elif masstype == 'm':
            self.density_delta = delta * cosmo.critical_density(redshift) * cosmo.Om(redshift)
        else:
            raise ValueError("Invalid masstype '" + str(masstype) +
                             "'. Should be one of 'virial','c', or 'm'")
        return self.density_delta

    @staticmethod
    def A_NFW(y):
        r"""
        Dimensionless volume integral of the NFW profile, used as an intermediate step in some
        calculations for this model.

        Notes
        -----

        Model formula:

        .. math:: A_{NFW} = [\ln(1+y) - \frac{y}{1+y}]
        """
        return np.log(1.0 + y) - (y / (1.0 + y))

    def _density_s(self, mass, concentration):
        """
        Calculate scale density of the NFW profile.
        """
        # Enforce default units
        if not isinstance(mass, u.Quantity):
            in_mass = u.Quantity(mass, u.M_sun)
        else:
            in_mass = mass

        # Calculate scale density
        # M_{200} = 4\pi \rho_{s} R_{s}^3 \left[\ln(1+c) - \frac{c}{1+c}\right].
        self.density_s = in_mass / (4.0 * np.pi * self._radius_s(in_mass, concentration) ** 3 *
                                    self.A_NFW(concentration))

        return self.density_s

    @property
    def rho_scale(self):
        r"""
        Scale density of the NFW profile. Often written in the literature as :math:`\rho_s`
        """
        return self.density_s

    def _radius_s(self, mass, concentration):
        """
        Calculate scale radius of the NFW profile.
        """
        # Enforce default units
        if not isinstance(mass, u.Quantity):
            in_mass = u.Quantity(mass, u.M_sun)
        else:
            in_mass = mass

        # Delta Mass is related to delta radius by
        # M_{200}=\frac{4}{3}\pi r_{200}^3 200 \rho_{c}
        # And delta radius is related to the NFW scale radius by
        # c = R / r_{\\rm s}
        self.radius_s = (((3.0 * in_mass) / (4.0 * np.pi * self.density_delta)) ** (
                          1.0 / 3.0)) / concentration

        # Set radial units to kiloparsec by default (unit will be rescaled by units of radius
        # in evaluate)
        return self.radius_s.to(u.kpc)

    @property
    def r_s(self):
        """
        Scale radius of the NFW profile.
        """
        return self.radius_s

    @property
    def r_virial(self):
        """
        Mass factor defined virial radius of the NFW profile (R200c for M200c, Rvir for Mvir, etc.).
        """
        return self.r_s * self.concentration

    @property
    def r_max(self):
        """
        Radius of maximum circular velocity.
        """
        return self.r_s * 2.16258

    @property
    def v_max(self):
        """
        Maximum circular velocity.
        """
        return self.circular_velocity(self.r_max)

    def circular_velocity(self, r):
        r"""
        Circular velocities of the NFW profile.

        Parameters
        ----------
        r : float or :class:`~astropy.units.Quantity`
            Radial position of velocity to be calculated for the NFW profile.

        Returns
        -------
        velocity : float or :class:`~astropy.units.Quantity`
            NFW profile circular velocity at location ``r``. The velocity units are:
            [km / s]

        Notes
        -----

        Model formula:

        .. math:: v_{circ}(r)^2 = \frac{1}{x}\frac{\ln(1+cx)-(cx)/(1+cx)}{\ln(1+c)-c/(1+c)}

        .. math:: x = r/r_s

        .. warning::

            Output values might contain ``nan`` and ``inf``.
        """
        # Enforce default units (if parameters are without units)
        if hasattr(r, "unit"):
            in_r = r
        else:
            in_r = u.Quantity(r, u.kpc)

        # Mass factor defined velocity (i.e. V200c for M200c, Rvir for Mvir)
        v_profile = np.sqrt(self.mass * const.G.to(in_r.unit**3 / (self.mass.unit * u.s**2)) /
                            self.r_virial)

        # Define reduced radius (r / r_{\\rm s})
        reduced_radius = in_r / self.r_virial.to(in_r.unit)

        # Circular velocity given by:
        # v^2=\frac{1}{x}\frac{\ln(1+cx)-(cx)/(1+cx)}{\ln(1+c)-c/(1+c)}
        # where x=r/r_{200}
        velocity = np.sqrt((v_profile**2 * self.A_NFW(self.concentration * reduced_radius)) /
                           (reduced_radius * self.A_NFW(self.concentration)))

        return velocity.to(u.km / u.s)

    @property
    def input_units(self):
        # The units for the 'r' variable should be a length (default kpc)
        return {self.inputs[0]: u.kpc}

    @property
    def return_units(self):
        # The units for the 'density' variable should be a matter density (default M_sun / kpc^3)
        if (self.mass.unit is None) and (self.input_units[self.inputs[0]] is None):
            return {self.outputs[0]: u.M_sun / u.kpc ** 3}
        elif (self.mass.unit is None):
            return {self.outputs[0]: u.M_sun / self.input_units[self.inputs[0]] ** 3}
        elif (self.input_units[self.inputs[0]] is None):
            return {self.outputs[0]: self.mass.unit / u.kpc ** 3}
        else:
            return {self.outputs[0]: self.mass.unit / self.input_units[self.inputs[0]] ** 3}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'mass': u.M_sun,
                "concentration": None,
                "redshift": None}
