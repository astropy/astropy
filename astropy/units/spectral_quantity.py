import astropy.units as u

__all__ = ['SpectralQuantity']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralQuantity.*']


SPECTRAL_UNITS = (u.Hz, u.m, u.J, (1 / u.m).unit, u.m / u.s)

DOPPLER_CONVENTIONS = {
    'radio': u.doppler_radio,
    'optical': u.doppler_optical,
    'relativistic': u.doppler_relativistic
}


class SpectralQuantity(u.SpecificTypeQuantity):
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

        kwargs['subok'] = True

        obj = super().__new__(cls, value, unit=unit, **kwargs)

        # The quantity machinery will drop the unit because type(value) !=
        #  SpectralCoord when passing in a Quantity object. Reassign the unit
        #  here to avoid this.
        # if isinstance(value, u.Quantity) and unit is None:
        #     obj._unit = value.unit

        # If we're initializing from an existing SpectralCoord, keep any
        # parameters that aren't being overridden
        if isinstance(value, SpectralQuantity):
            if doppler_rest is None:
                doppler_rest = value.doppler_rest
            if doppler_convention is None:
                doppler_convention = value.doppler_convention

        obj._doppler_rest = doppler_rest
        obj._doppler_convention = doppler_convention

        return obj

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)
        self._doppler_rest = getattr(obj, '_doppler_rest', None)
        self._doppler_convention = getattr(obj, '_doppler_convention', None)

    def __quantity_subclass__(self, unit):
        """
        Overridden by subclasses to change what kind of view is
        created based on the output unit of an operation.
        """
        return SpectralQuantity, True

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
    @u.quantity_input(value=SPECTRAL_UNITS)
    def doppler_rest(self, value):
        """
        New rest value needed for velocity-space conversions.

        Parameters
        ----------
        value : `~astropy.units.Quantity`
            Rest value.
        """
        if self._doppler_rest is None:
            self._doppler_rest = value
        else:
            raise ValueError("doppler_rest has already been set, and cannot "
                             "be changed. Use the ``to`` method to convert "
                             "the spectral values(s) to use a different "
                             "rest value")

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
        if value is not None and value not in DOPPLER_CONVENTIONS:
            raise ValueError("doppler_convention should be one of {0}".format('/'.join(sorted(DOPPLER_CONVENTIONS))))

        if self._doppler_convention is None:
            self._doppler_convention = value
        else:
            raise ValueError("doppler_convention has already been set, and cannot "
                             "be changed. Use the ``to`` method to convert "
                             "the spectral values(s) to use a different "
                             "convention")

    @u.quantity_input(doppler_rest=SPECTRAL_UNITS)
    def to(self, unit,
           equivalencies=[],
           doppler_rest=None,
           doppler_convention=None):
        """
        Overloaded parent ``to`` method to provide parameters for defining
        rest value and pre-defined conventions for unit transformations.

        Parameters
        ----------
        doppler_rest : `~astropy.units.Quantity`, optional
            The rest value used in the velocity space conversions. Providing
            the value here will set the value stored on the `SpectralCoord`
            instance.
        doppler_convention : {'relativistic', 'optical', 'radio'}, optional
            The velocity convention to use during conversions. Providing the
            value here will set the value stored on the `SpectralCoord`
            instance.

        Returns
        -------
        `SpectralCoord`
            New spectral coordinate object with data converted to the new unit.
        """

        # FIXME: need to consider case where doppler equivalency is passed in
        # equivalencies list, or is u.spectral equivalency is already passed

        if doppler_rest is None:
            doppler_rest = self._doppler_rest

        if doppler_convention is None:
            doppler_convention = self._doppler_convention
        elif doppler_convention not in DOPPLER_CONVENTIONS:
            raise ValueError("doppler_convention should be one of {0}".format('/'.join(sorted(DOPPLER_CONVENTIONS))))

        if self.unit.is_equivalent(u.km / u.s) and unit.is_equivalent(u.km / u.s):

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

            freq = super().to(u.Hz, equivalencies=equivalencies + u.spectral() + vel_equiv1)

            vel_equiv2 = DOPPLER_CONVENTIONS[doppler_convention](doppler_rest)

            result = freq.to(unit, equivalencies=equivalencies + u.spectral() + vel_equiv2)

        else:

            additional_equivalencies = u.spectral()

            if self.unit.is_equivalent(u.km / u.s) or unit.is_equivalent(u.km / u.s):

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
