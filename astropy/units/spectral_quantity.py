import numpy as np
from . import si
from . import equivalencies as eq
from .quantity import SpecificTypeQuantity
from .decorators import quantity_input
from astropy.constants import c

__all__ = ['SpectralQuantity']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralQuantity.*']


SPECTRAL_UNITS = (si.Hz, si.m, si.J, si.m ** -1, si.m / si.s)

DOPPLER_CONVENTIONS = {
    'radio': eq.doppler_radio,
    'optical': eq.doppler_optical,
    'relativistic': eq.doppler_relativistic
}


class SpectralQuantity(SpecificTypeQuantity):
    """
    One or more value(s) with spectral units.

    The spectral units should be those for frequencies, wavelengths, energies,
    wavenumbers, or velocities (interpreted as Doppler velocities relative to a
    rest spectral value). The advantage of using this class over the regular
    `~astropy.units.Quantity` class is that in `SpectralQuantity`, the
    ``u.spectral`` equivalency is enabled by default (allowing automatic
    conversion between spectral units), and a preferred Doppler rest value and
    convention can be stored for easy conversion to/from velocities.

    Parameters
    ----------
    value : ndarray or `~astropy.units.Quantity` or `SpectralQuantity`
        Spectral axis data values.
    unit : str or `~astropy.units.Unit`
        Unit for the given data.
    doppler_rest : `~astropy.units.Quantity`, optional
        The rest value to use for conversions from/to velocities
    doppler_convention : str, optional
        The convention to use when converting the spectral data to/from
        velocities.
    """

    _equivalent_unit = SPECTRAL_UNITS

    _include_easy_conversion_members = True

    def __new__(cls, value, unit=None,
                doppler_rest=None, doppler_convention=None,
                **kwargs):

        obj = super().__new__(cls, value, unit=unit, **kwargs)

        # If we're initializing from an existing SpectralQuantity, keep any
        # parameters that aren't being overridden
        if doppler_rest is None:
            doppler_rest = getattr(value, 'doppler_rest', None)
        if doppler_convention is None:
            doppler_convention = getattr(value, 'doppler_convention', None)

        obj._doppler_rest = doppler_rest
        obj._doppler_convention = doppler_convention

        return obj

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)
        self._doppler_rest = getattr(obj, '_doppler_rest', None)
        self._doppler_convention = getattr(obj, '_doppler_convention', None)

    @property
    def doppler_rest(self):
        """
        The rest value of the spectrum used for transformations to/from
        velocity space.

        Returns
        -------
        `~astropy.units.Quantity`
            Rest value as an astropy `~astropy.units.Quantity` object.
        """
        return self._doppler_rest

    @doppler_rest.setter
    @quantity_input(value=SPECTRAL_UNITS)
    def doppler_rest(self, value):
        """
        New rest value needed for velocity-space conversions.

        Parameters
        ----------
        value : `~astropy.units.Quantity`
            Rest value.
        """
        if self._doppler_rest is not None:
            raise AttributeError("doppler_rest has already been set, and cannot "
                                 "be changed. Use the ``to`` method to convert "
                                 "the spectral values(s) to use a different "
                                 "rest value")
        self._doppler_rest = value

    @property
    def doppler_convention(self):
        """
        The defined convention for conversions to/from velocity space.

        Returns
        -------
        str
            One of 'optical', 'radio', or 'relativistic' representing the
            equivalency used in the unit conversions.
        """
        return self._doppler_convention

    @doppler_convention.setter
    def doppler_convention(self, value):
        """
        New velocity convention used for velocity space conversions.

        Parameters
        ----------
        value

        Notes
        -----
        More information on the equations dictating the transformations can be
        found in the astropy documentation [1]_.

        References
        ----------
        .. [1] Astropy documentation: https://docs.astropy.org/en/stable/units/equivalencies.html#spectral-doppler-equivalencies

        """

        if self._doppler_convention is not None:
            raise AttributeError("doppler_convention has already been set, and cannot "
                                 "be changed. Use the ``to`` method to convert "
                                 "the spectral values(s) to use a different "
                                 "convention")

        if value is not None and value not in DOPPLER_CONVENTIONS:
            raise ValueError("doppler_convention should be one of {0}".format('/'.join(sorted(DOPPLER_CONVENTIONS))))

        self._doppler_convention = value

    @quantity_input(doppler_rest=SPECTRAL_UNITS)
    def to(self, unit,
           equivalencies=[],
           doppler_rest=None,
           doppler_convention=None):
        """
        Return a new `~astropy.units.SpectralQuantity` object with the specified unit.

        By default, the ``spectral`` equivalency will be enabled, as well as
        one of the Doppler equivalencies if converting to/from velocities.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to convert to. Must be
            an `~astropy.units.UnitBase` object or a string parseable
            by the `~astropy.units` package, and should be a spectral unit.
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            If not provided or ``[]``, class default equivalencies will be used
            (none for `~astropy.units.Quantity`, but may be set for subclasses)
            If `None`, no equivalencies will be applied at all, not even any
            set globally or within a context.
        doppler_rest : `~astropy.units.Quantity`, optional
            The rest value used when converting to/from velocities. This will
            also be set at an attribute on the output
            `~astropy.units.SpectralQuantity`.
        doppler_convention : {'relativistic', 'optical', 'radio'}, optional
            The Doppler convention used when converting to/from velocities.
            This will also be set at an attribute on the output
            `~astropy.units.SpectralQuantity`.

        Returns
        -------
        `SpectralQuantity`
            New spectral coordinate object with data converted to the new unit.
        """

        # If equivalencies is explicitly set to None, we should just use the
        # default Quantity.to with equivalencies also set to None
        if equivalencies is None:
            return super().to(unit, equivalencies=None)

        # FIXME: need to consider case where doppler equivalency is passed in
        # equivalencies list, or is u.spectral equivalency is already passed

        if doppler_rest is None:
            doppler_rest = self._doppler_rest

        if doppler_convention is None:
            doppler_convention = self._doppler_convention
        elif doppler_convention not in DOPPLER_CONVENTIONS:
            raise ValueError("doppler_convention should be one of {0}".format('/'.join(sorted(DOPPLER_CONVENTIONS))))

        if self.unit.is_equivalent(si.km / si.s) and unit.is_equivalent(si.km / si.s):

            # Special case: if the current and final units are both velocity,
            # and either the rest value or the convention are different, we
            # need to convert back to frequency temporarily.

            if doppler_convention is not None and self._doppler_convention is None:
                raise ValueError("Original doppler_convention not set")

            if doppler_rest is not None and self._doppler_rest is None:
                raise ValueError("Original doppler_rest not set")

            if doppler_rest is None and doppler_convention is None:
                return super().to(unit, equivalencies=equivalencies)
            elif (doppler_rest is None) is not (doppler_convention is None):
                raise ValueError("Either both or neither doppler_rest and "
                                 "doppler_convention should be defined for "
                                 "velocity conversions")

            vel_equiv1 = DOPPLER_CONVENTIONS[self._doppler_convention](self._doppler_rest)

            freq = super().to(si.Hz, equivalencies=equivalencies + vel_equiv1)

            vel_equiv2 = DOPPLER_CONVENTIONS[doppler_convention](doppler_rest)

            result = freq.to(unit, equivalencies=equivalencies + vel_equiv2)

        else:

            additional_equivalencies = eq.spectral()

            if self.unit.is_equivalent(si.km / si.s) or unit.is_equivalent(si.km / si.s):

                if doppler_convention is None:
                    raise ValueError("doppler_convention not set, cannot convert to/from velocities")

                if doppler_rest is None:
                    raise ValueError("doppler_rest not set, cannot convert to/from velocities")

                additional_equivalencies += DOPPLER_CONVENTIONS[doppler_convention](doppler_rest)

            result = super().to(unit, equivalencies=equivalencies + additional_equivalencies)

        result._doppler_rest = doppler_rest
        result._doppler_convention = doppler_convention

        return result

    def to_value(self, *args, **kwargs):
        return self.to(*args, **kwargs).value

    def _apply_relativistic_doppler_shift(self, velocity):
        """
        Given a `SpectralQuantity` and a velocity, return a new `SpectralQuantity`
        that is Doppler shifted by this amount.

        Note that the Doppler shift applied is the full relativistic one, so
        `SpectralQuantity` currently expressed in velocity and not using the
        relativistic convention will temporarily be converted to use the
        relativistic convention while the shift is applied.

        Positive velocities are assumed to redshift the spectral quantity,
        while negative velocities blueshift the spectral quantity.
        """

        # NOTE: we deliberately don't keep sub-classes of SpectralQuantity intact
        # since we can't guarantee that their metadata would be correct/consistent.
        squantity = self.view(SpectralQuantity)

        beta = velocity / c
        doppler_factor = np.sqrt((1 + beta) / (1 - beta))

        if squantity.unit.is_equivalent(si.m):  # wavelength
            return squantity * doppler_factor
        elif squantity.unit.is_equivalent(si.Hz) or squantity.unit.is_equivalent(si.eV) or squantity.unit.is_equivalent(1 / si.m):
            return squantity / doppler_factor
        elif squantity.unit.is_equivalent(si.km / si.s):  # velocity
            if squantity.doppler_convention is None:
                raise ValueError('doppler_convention is not set, so unsure how to apply Doppler shift')
            return (squantity.to(si.Hz) / doppler_factor).to(squantity.unit)
        else:
            raise RuntimeError(f"Unexpected units in velocity shift: {squantity.unit}. "
                               "This should not happen, so please report this in the "
                               "astropy issue tracker!")
