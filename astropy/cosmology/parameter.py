# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
from astropy.utils.decorators import classproperty

__all__ = ["Parameter"]


class Parameter:
    r"""Cosmological parameter (descriptor).

    Should only be used with a :class:`~astropy.cosmology.Cosmology` subclass.

    Parameters
    ----------
    fvalidate : callable[[object, object, Any], Any] or str, optional
        Function to validate the Parameter value from instances of the
        cosmology class. If "default", uses default validator to assign units
        (with equivalencies), if Parameter has units.
        For other valid string options, see ``Parameter._registry_validators``.
        'fvalidate' can also be set through a decorator with
        :meth:`~astropy.cosmology.Parameter.validator`.
    doc : str or None, optional
        Parameter description.
    unit : unit-like or None (optional, keyword-only)
        The `~astropy.units.Unit` for the Parameter. If None (default) no
        unit as assumed.
    equivalencies : `~astropy.units.Equivalency` or sequence thereof
        Unit equivalencies for this Parameter.
    fmt : str (optional, keyword-only)
        `format` specification, used when making string representation
        of the containing Cosmology.
        See https://docs.python.org/3/library/string.html#formatspec
    derived : bool (optional, keyword-only)
        Whether the Parameter is 'derived', default `False`.
        Derived parameters behave similarly to normal parameters, but are not
        sorted by the |Cosmology| signature (probably not there) and are not
        included in all methods. For reference, see ``Ode0`` in
        ``FlatFLRWMixin``, which removes :math:`\Omega_{de,0}`` as an
        independent parameter (:math:`\Omega_{de,0} \equiv 1 - \Omega_{tot}`).

    Examples
    --------
    For worked examples see :class:`~astropy.cosmology.FLRW`.
    """

    _registry_validators = {}

    def __init__(self, fvalidate="default", doc=None, *,
                 unit=None, equivalencies=[], fmt=".3g", derived=False):
        # parse registered fvalidate
        if callable(fvalidate):
            pass
        elif fvalidate in self._registry_validators:
            fvalidate = self._registry_validators[fvalidate]
        elif isinstance(fvalidate, str):
            raise ValueError("`fvalidate`, if str, must be in "
                             f"{self._registry_validators.keys()}")
        else:
            raise TypeError("`fvalidate` must be a function or "
                            f"{self._registry_validators.keys()}")

        self.__doc__ = doc
        self._fvalidate = fvalidate

        # units stuff
        self._unit = u.Unit(unit) if unit is not None else None
        self._equivalencies = equivalencies

        # misc
        self._fmt = str(fmt)
        self._derived = derived

    def __set_name__(self, cosmo_cls, name):
        # attribute name
        self._attr_name = name
        self._attr_name_private = "_" + name

    @property
    def name(self):
        """Parameter name."""
        return self._attr_name

    @property
    def unit(self):
        """Parameter unit."""
        return self._unit

    @property
    def equivalencies(self):
        """Equivalencies used when initializing Parameter."""
        return self._equivalencies

    @property
    def format_spec(self):
        """String format specification."""
        return self._fmt

    @property
    def derived(self):
        """Whether the Parameter is derived; true parameters are not."""
        return self._derived

    # -------------------------------------------
    # descriptor and property-like methods

    def __get__(self, cosmology, cosmo_cls=None):
        # get from class
        if cosmology is None:
            return self
        return getattr(cosmology, self._attr_name_private)

    def __set__(self, cosmology, value):
        """Allows attribute setting once. Raises AttributeError subsequently."""
        # raise error if setting 2nd time.
        if hasattr(cosmology, self._attr_name_private):
            raise AttributeError("can't set attribute")

        # validate value, generally setting units if present
        value = self.validate(cosmology, value)
        setattr(cosmology, self._attr_name_private, value)

    # -------------------------------------------
    # validate value

    @property
    def fvalidate(self):
        """Function to validate a potential value of this Parameter.."""
        return self._fvalidate

    def validator(self, fvalidate):
        """Make new Parameter with custom ``fvalidate``.

        Note: ``Parameter.fvalidator`` must be the top-most descriptor decorator.

        Parameters
        ----------
        fvalidate : callable[[type, type, Any], Any]

        Returns
        -------
        `~astropy.cosmology.Parameter`
            Copy of this Parameter but with custom ``fvalidate``.
        """
        desc = type(self)(fvalidate=fvalidate,
                          doc=self.__doc__, fmt=self.format_spec,
                          unit=self.unit, equivalencies=self.equivalencies,
                          derived=self.derived)
        return desc

    def validate(self, cosmology, value):
        """Run the validator on this Parameter.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.Cosmology` instance
        value : Any
            The object to validate.

        Returns
        -------
        Any
            The output of calling ``fvalidate(cosmology, self, value)``
            (yes, that parameter order).
        """
        return self.fvalidate(cosmology, self, value)

    @classmethod
    def register_validator(cls, key, fvalidate=None):
        """Decorator to register a new kind of validator function.

        Parameters
        ----------
        key : str
        fvalidate : callable[[object, object, Any], Any] or None, optional
            Value validation function.

        Returns
        -------
        ``validator`` or callable[``validator``]
            if validator is None returns a function that takes and registers a
            validator. This allows ``register_validator`` to be used as a
            decorator.
        """
        if key in cls._registry_validators:
            raise KeyError(f"validator {key!r} already registered with Parameter.")

        # fvalidate directly passed
        if fvalidate is not None:
            cls._registry_validators[key] = fvalidate
            return fvalidate

        # for use as a decorator
        def register(fvalidate):
            """Register validator function.

            Parameters
            ----------
            fvalidate : callable[[object, object, Any], Any]
                Validation function.

            Returns
            -------
            ``validator``
            """
            cls._registry_validators[key] = fvalidate
            return fvalidate

        return register

    # -------------------------------------------

    def __repr__(self):
        return f"<Parameter {self._attr_name!r} at {hex(id(self))}>"


# ===================================================================
# Built-in validators


@Parameter.register_validator("default")
def _validate_with_unit(cosmology, param, value):
    """
    Default Parameter value validator.
    Adds/converts units if Parameter has a unit.
    """
    if param.unit is not None:
        with u.add_enabled_equivalencies(param.equivalencies):
            value = u.Quantity(value, param.unit)
    return value


@Parameter.register_validator("float")
def _validate_to_float(cosmology, param, value):
    """Parameter value validator with units, and converted to float."""
    value = _validate_with_unit(cosmology, param, value)
    return float(value)


@Parameter.register_validator("scalar")
def _validate_to_scalar(cosmology, param, value):
    """"""
    value = _validate_with_unit(cosmology, param, value)
    if not value.isscalar:
        raise ValueError(f"{param.name} is a non-scalar quantity")
    return value


@Parameter.register_validator("non-negative")
def _validate_non_negative(cosmology, param, value):
    """Parameter value validator where value is a positive float."""
    value = _validate_to_float(cosmology, param, value)
    if value < 0.0:
        raise ValueError(f"{param.name} cannot be negative.")
    return value
